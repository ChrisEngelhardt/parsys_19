#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>

int SIZE = 100;

int** a;
int** b;
int** serialResult;
int** mpResult;

void genMat(){
	a=malloc(sizeof(int*)*SIZE);
	b=malloc(sizeof(int*)*SIZE);
	serialResult=malloc(sizeof(int*)*SIZE);
	mpResult=malloc(sizeof(int*)*SIZE);
	for(int i=0;i<SIZE;i++){
		a[i]=malloc(sizeof(int)*SIZE);
		b[i]=malloc(sizeof(int)*SIZE);
		serialResult[i]=malloc(sizeof(int)*SIZE);
		mpResult[i]=malloc(sizeof(int)*SIZE);
		for(int j=0;j<SIZE;j++){
			a[i][j]=rand()%100;
			b[i][j]=rand()%100;
			serialResult[i][j]=0;
			mpResult[i][j]=0;
		}
	}
}

void matMulSerial(){
	for(int i=0;i<SIZE;i++){
		for(int j=0;j<SIZE;j++){
			for(int k=0;k<SIZE;k++){
				serialResult[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
}

void verify(){
	for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++){
            assert(serialResult[i][j]==mpResult[i][j]);
		}
	}
}

void matMulMP(){
	#pragma omp parallel for shared(a,b,mpResult)
	for(int i=0;i<SIZE;i++){
		for(int j=0;j<SIZE;j++){
			for(int k=0;k<SIZE;k++){
				mpResult[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
}

void freeMat(){
	for(int i=0;i<SIZE;i++){
		free(a[i]);
		free(b[i]);
		free(serialResult[i]);
		free(mpResult[i]);
	}
	free(a);
	free(b);
	free(serialResult);
	free(mpResult);
}


int main (int argc, char* argv[]){

    if (argc > 1) SIZE = atoi(argv[1]);

	int THREADCOUNT = omp_get_max_threads();
    genMat();
    matMulSerial();

    double startTime = omp_get_wtime();
    matMulMP();
	printf("%d; ", THREADCOUNT);
	printf("%d; ", SIZE);
  	printf("%lf\n", omp_get_wtime() - startTime);
    verify();

	freeMat();
}