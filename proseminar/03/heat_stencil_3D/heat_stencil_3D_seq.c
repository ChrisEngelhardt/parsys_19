#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double value_t;

#define RESOLUTION 50

// -- vector utilities --

typedef value_t *Vector;

typedef value_t **Matrix;

typedef value_t ***Cube;


Vector createVector(int N);

Cube createCube(int x, int y, int z);

void releaseVector(Vector m);

void releaseCube(Cube m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
  //Using for time mesurement
  MPI_Init(0,0);
  
  double start = MPI_Wtime();
  
  // 'parsing' optional input parameter = problem size
  int N = 50;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 10;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Cube A = createCube(N,N,N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++){
      for(int k = 0; k < N; k++){
        A[i][j][k] = 273; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = N / 4;
  int source_z = N / 4;

  
  A[source_x][source_y][source_z] = 273 + 60;
  
  // ---------- compute ----------

  // create a second buffer for the computation
  Cube B = createCube(N,N,N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long i = 0; i < N; i++) {
      for(long long j = 0; j < N; j++){
        for(long long k = 0; k < N; k++){

          // center stays constant (the heat is still on)
          if (i == source_x && j == source_y && k == source_z) {
            B[i][j][k] = A[i][j][k];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i][j][k];

          // get temperatures of adjacent cells
          value_t tl = (j != 0) ? A[i][j - 1][k] : tc;
          value_t tr = (j != N - 1) ? A[i][j + 1][k] : tc;
          value_t tb = (i != 0) ? A[i-1][j][k]: tc;
          value_t tt = (i != N - 1) ? A[i+1][j][k] : tc;
          value_t tba = (k != 0) ? A[i][j][k-1]: tc;
          value_t tfr = (k != N - 1) ? A[i][j][k+1] : tc;


          // compute new temperature at current position
          B[i][j][k] = tc + 0.15 * (tl + tr + tb + tt + tba + tfr + (-6 * tc));
        }
      }
    }

    // swap matrices (just pointers, not content)
    Cube H = A;
    A = B;
    B = H;
  }

  releaseCube(B,N);
  
  int success = 1;
  for (long long i = 0; i < N; i++) {
    for(long long j = 0; j < N; j++){
      for(long long k = 0; j < N; j++){
        value_t temp = A[i][j][k];
        if (273 <= temp && temp <= 273 + 60)
          continue;
        success = 0;
        break;
      }
    }
  }
  
  printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  
  printf("%lf",MPI_Wtime() - start);

  
  // ---------- cleanup ----------

  releaseCube(A,N);
  
  MPI_Finalize();

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

Cube createCube(int x, int y, int z) {
  Cube cube = malloc(sizeof(Matrix) * x);
  
  for(int i = 0; i < x; i ++){
    cube[i] = malloc(sizeof(Vector) * y);
  }
  
  cube[0][0] = malloc(sizeof(value_t) * x * y * z);

  for(int i = 0; i < x; i++){
    for(int j = 0; j < y; j++){
      cube[i][j] = &cube[0][0][i*y*z+j*z];
    }
  }
  
  return cube;
}
void releaseVector(Vector m) { free(m); }

void releaseCube(Cube c, int x) {
  free(c[0][0]);
  
  for(int i = 0; i < x; i ++){
    free(c[i]);
  }
  free(c);
}