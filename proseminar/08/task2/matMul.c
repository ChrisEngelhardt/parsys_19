#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>

#define SIZE 2000

int a[SIZE][SIZE];
int b[SIZE][SIZE];
int serialResult[SIZE][SIZE];
int mpResult[SIZE][SIZE];
int THREADCOUNT = 4;

void genMat(){
	int i,j;
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			a[i][j]=rand()%100;
			b[i][j]=rand()%100;
		}
	}
}

void matMulSerial(){
	int i,j,k;
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			for(k=0;k<SIZE;k++){
				serialResult[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
}

void verify(){
	int i,j;
	for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
            assert(serialResult[i][j]==mpResult[i][j]);
}

void matMulMP(){
	int i,j,k;
	#pragma omp parallel for private(i,j,k) shared(a,b,mpResult) num_threads(THREADCOUNT) 
	for(i=0;i<SIZE;i++)
        for( j=0;j<SIZE;j++)
            for(k=0;k<SIZE;k++)
                mpResult[i][j]+=a[i][k]*b[k][j];

}


int main (int argc, char* argv[]){

    if (argc > 1) THREADCOUNT = atoi(argv[1]);
    genMat();
    matMulSerial();

    double startTime = omp_get_wtime();
    matMulMP();
    printf("%lf; ", omp_get_wtime() - startTime);
    verify();
}