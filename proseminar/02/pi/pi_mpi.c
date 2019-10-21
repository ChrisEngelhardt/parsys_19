#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>


#define CIRCLELENGTH 1
#define TAG 42

int main (int argc, char* argv[]){
    int numProcs, myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (argc < 2){
        printf("Please provide number of samples\n");
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    srand((unsigned int)time(NULL)*myrank);
    long inside = 0;
    long iterations = atol(argv[1]);
    printf("%ld samples\n", iterations);
    
    
    
    for(long i = 0; i < iterations/numProcs; i++){
        double point[2];
        point[0] = (double)rand()/(RAND_MAX);   
        point[1] = (double)rand()/(RAND_MAX);   
        if (point[0]*point[0]+point[1]*point[1] < CIRCLELENGTH){
            inside ++;
        }
    }

    if(myrank == 0){
        for(int i = 0; i < numProcs -1 ;i++){
            long recvBuff = 0;
            MPI_Recv(&recvBuff, 1, MPI_LONG, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            inside += recvBuff;
        } 
        printf("Pi is approx. %.5f\n",(double) inside/((iterations/numProcs)*numProcs) * 4);
    }else{
        MPI_Send(&inside, 1, MPI_LONG, 0, TAG, MPI_COMM_WORLD);
    } 

    MPI_Finalize();
    return EXIT_SUCCESS;
}