#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define CIRCLELENGTH 1

int main (int argc, char* argv[]){

    if (argc < 2){
        printf("Please provide number of samples\n");
        return EXIT_SUCCESS;
    }

    srand((unsigned int)time(NULL));
    int inside = 0;
    long iterations = atol(argv[1]);
    printf("%ld samples\n", iterations);
    for(int i = 0; i < iterations; i++){
        double point[2];
        point[0] = (double)rand()/(RAND_MAX);   
        point[1] = (double)rand()/(RAND_MAX);   
        if (point[0]*point[0]+point[1]*point[1] < CIRCLELENGTH){
            inside ++;
        }
      
    }

    printf("Pi is approx. %.5f\n",(double) inside/iterations * 4);
    return EXIT_SUCCESS;
}