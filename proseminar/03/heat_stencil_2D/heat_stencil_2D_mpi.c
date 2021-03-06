#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double value_t;

#define TAG_L 42
#define TAG_R 43
#define TAG_A 44


#define RESOLUTION 50

// -- vector utilities --

typedef value_t *Vector;

typedef value_t **Matrix;

Vector createVector(int N);

Matrix createMatrix(int x, int y);

void releaseVector(Vector m);

void releaseMatrix(Matrix m);

void printTemperature(Vector m, int N);

void printTemperatureMatrix(Matrix m, int x, int y);

// -- simulation code ---

int main(int argc, char **argv) {
  int numProcs, myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  double start = MPI_Wtime();
  
  int lastProcess = numProcs-1;
  
  
  // 'parsing' optional input parameter = problem size
  int N = 200;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 10;
  if(myrank == lastProcess){
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Matrix A = createMatrix(N,N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++){
      A[i][j] = 273; // temperature is 0° C everywhere (273 K)
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = N / 4;
  
  A[source_x][source_y] = 273 + 60;

  if(myrank == lastProcess){
    printf("Initial:\n");
    printTemperatureMatrix(A, N, N);
    printf("\n");
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Matrix B = createMatrix(N,N);
  
  long long NPerSlot = N/numProcs;

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long i = NPerSlot*myrank; i < NPerSlot*(myrank+1) || (myrank==lastProcess && i < N); i++) {
      for(long long j = 0; j < N; j++){
        // center stays constant (the heat is still on)
        if (i == source_x && j == source_y) {
          B[i][j] = A[i][j];
          continue;
        }

        // get temperature at current position
        value_t tc = A[i][j];

        // get temperatures of adjacent cells
        value_t tl = (j != 0) ? A[i][j - 1] : tc;
        value_t tr = (j != N - 1) ? A[i][j + 1] : tc;
        value_t tb = (i != 0) ? A[i-1][j] : tc;
        value_t tt = (i != N - 1) ? A[i+1][j] : tc;

        // compute new temperature at current position
        B[i][j] = tc + 0.2 * (tl + tr + tb + tt + (-4 * tc));
      }
    }

    // swap matrices (just pointers, not content)
    Matrix H = A;
    A = B;
    B = H;
    
    MPI_Request requestL;
    MPI_Request requestR;
    //Send
    if(myrank != 0){
      MPI_Isend(&A[NPerSlot*myrank][0], N, MPI_DOUBLE, myrank-1, TAG_L, MPI_COMM_WORLD, &requestL);
    }
    
    if(myrank != lastProcess){
      MPI_Isend(&A[NPerSlot*(myrank+1)-1][0], N, MPI_DOUBLE, myrank+1, TAG_R, MPI_COMM_WORLD, &requestR);
    }
    
    //Recive
    if(myrank != 0){
      MPI_Recv(&A[NPerSlot*myrank-1][0], N, MPI_DOUBLE, myrank-1, TAG_R, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(myrank != lastProcess){
      MPI_Recv(&A[NPerSlot*(myrank+1)][0], N, MPI_DOUBLE, myrank+1, TAG_L, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    
    //Waiting for completion
    if(myrank != 0){
      MPI_Wait(&requestL, MPI_STATUS_IGNORE);
    }
    if(myrank != lastProcess){
      MPI_Wait(&requestR, MPI_STATUS_IGNORE);
    }
  }

  releaseMatrix(B);
  
  MPI_Gather(A[NPerSlot*myrank], NPerSlot*N, MPI_DOUBLE, A[0], NPerSlot*N, MPI_DOUBLE, lastProcess, MPI_COMM_WORLD);

  // ---------- check ----------

  if(myrank == lastProcess){
    printf("Final:\n");
    printTemperatureMatrix(A, N, N);
    printf("\n");
  }
  int success = 1;
  for (long long i = 0; i < N; i++) {
    for(long long j = 0; j < N; j++){
      value_t temp = A[i][j];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }
  }

  if(myrank == lastProcess){
    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  }

  // ---------- cleanup ----------

  releaseMatrix(A);
  
  if(myrank == lastProcess){
    printf("%lf",MPI_Wtime() - start);
  }
  
  MPI_Finalize();

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

Matrix createMatrix(int x, int y) {
  Matrix matrix = malloc(sizeof(Vector) * x);
  matrix[0] = malloc(sizeof(value_t) * x * y);

  for(int i = 1; i < x; i++){
    matrix[i] = &matrix[0][i*y];
  }
  
  return matrix;
}
void releaseVector(Vector m) { free(m); }

void releaseMatrix(Matrix m) {
  free(m[0]);
  free(m);
}


void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}

void printTemperatureMatrix(Matrix m, int x, int y){
  printf("\e[1;1H\e[2J");
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW_x = x / W;
  int sW_y = y / W;

  // room

  // actual room
  for (int i = 0; i < W; i++) {
    // left wall
    printf("X");
    for (int j = 0; j < W; j++){
      value_t max_t = 0;
      // get max temperature in this tile
      for (int x = sW_x * i; x < sW_x * i + sW_x; x++) {
        for (int y = sW_y * j; y < sW_y * j + sW_y; y++) {
          max_t = (max_t < m[x][y]) ? m[x][y] : max_t;
        }
      }
      
      value_t temp = max_t;

      // pick the 'color'
      int c = ((temp - min) / (max - min)) * numColors;
      c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

      // print the average temperature
      printf("%c", colors[c]);
    }
    
    // right wall
    printf("X\n");
  }

}
