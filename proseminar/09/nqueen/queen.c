#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int solutions = 0;
int n = 8;

void setQueen(int positions[], int row, int col) {
	for(int i=0; i<row; i++) {
		// Check if position is safe
		if ((positions[i]==col) || (abs(positions[i]-col) == (row-i))) {
			return;
		}
	}
	// col is ok, set queen.
	positions[row]=col;
	if(row==n-1) {
		solutions++;
	} else {
		for(int i=0; i<n; i++) {
			setQueen(positions, row+1, i);
		}
	}
}

int main(int argc, char* argv[]) {
	n = (argc > 1) ? atoi(argv[1]) : 8;
	int* array = malloc(sizeof(int)*n);
	
	double start_time, end_time;
	
	start_time = omp_get_wtime();
	
	#pragma omp parallel for
	for(int i=0; i<n; i++)
		setQueen(array, 0, i);
	
	end_time = omp_get_wtime();
		
	free(array);
	
	printf("The execution time is %g sec\n", end_time - start_time);
	printf("Number of found solutions is %d\n", solutions);
		
	return 0;
}