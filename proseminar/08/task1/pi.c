#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <omp.h>


#define CIRCLELENGTH 1


int main (int argc, char* argv[]){
 	
    double startTime = omp_get_wtime();
    long iterations = 1000;
    int THREADCOUNT = 4;

    if (argc > 1) iterations = atoi(argv[1]);
    if (argc > 2) THREADCOUNT = atoi(argv[2]);

    srand((unsigned int)time(NULL));
    
    printf("%ld; ", iterations);

    long inside = 0;
    #pragma omp parallel for num_threads(THREADCOUNT) reduction(+:inside)
    for(long i = 0; i < iterations; i++){
        double point[2];
        point[0] = (double)rand()/(RAND_MAX);   
        point[1] = (double)rand()/(RAND_MAX);   
        if (point[0]*point[0]+point[1]*point[1] < CIRCLELENGTH){
            inside ++;
        }
    }
    printf("%lf; ", omp_get_wtime() - startTime);
    printf("%.5f\n",(double) inside/iterations * 4);
}