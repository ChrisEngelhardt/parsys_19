#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>
#define CIRCLELENGTH 1


int main (int argc, char* argv[]){
 	
    double startTime = omp_get_wtime();
    long iterations = 1000;
    int THREADCOUNT = omp_get_max_threads();

    if (argc > 1) iterations = atoi(argv[1]);

    
    //printf("%ld; ", iterations);

    long inside = 0;
    #pragma omp parallel
    {
        int r = omp_get_thread_num();

        #pragma omp for reduction(+:inside) 
        for(long i = 0; i < iterations; i++){
            double point[2];
            point[0] = (double)rand_r(&r)/(RAND_MAX);   
            point[1] = (double)rand_r(&r)/(RAND_MAX);   
            if (point[0]*point[0]+point[1]*point[1] < CIRCLELENGTH){
                inside ++;
            }
        }
    }
    printf("%d; ", THREADCOUNT);
    printf("%ld; ", iterations);
    printf("%lf\n", omp_get_wtime() - startTime);
    //printf("%.5f\n",(double) inside/iterations * 4);
}