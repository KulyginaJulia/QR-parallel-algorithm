#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>

void generator(double* MatrixA, double* MatrixB, double* MatrixC, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int tmp = rand() % 100;
			MatrixA[i*N + j] = (double)(tmp);
			tmp = rand() % 100;
			MatrixB[i*N + j] = (double)(tmp);
			MatrixC[i*N + j] = 0;
		}
	}
}
void multiplication_block(double* MatrixA, double* MatrixB, double* MatrixC, int Size) { 
	int BlockSize = 250;
	int GridSize = (int)(Size / (double)(BlockSize));
	#pragma omp parallel for
	for (int n = 0; n<GridSize; n++)
		for (int m = 0; m<GridSize; m++)
			for (int iter = 0; iter<GridSize; iter++)
				for (int i = n * BlockSize; i<(n + 1)*BlockSize; i++)
					for (int j = m * BlockSize; j<(m + 1)*BlockSize; j++)
						for (int k = iter * BlockSize; k<(iter + 1)*BlockSize; k++) 
							MatrixC[i*Size + j] += MatrixA[i*Size + k] * MatrixB[k*Size + j];
}

void print(double* MatrixA, double* MatrixB, double* MatrixC, int Size) {
	printf("Matrix A:\n");
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			printf("%.6f  ", MatrixA[i*Size + j]);
		}
		printf("\n");
	}

	printf("Matrix B:\n");
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			printf("%.6f  ", MatrixB[i*Size + j]);
		}
		printf("\n");
	}
	printf("Matrix C:\n");
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			printf("%.6f  ", MatrixC[i*Size + j]);
		}
		printf("\n");
	}
}
int main(int argc, char **argv) {

	int t = 0;
	int Size = 0; 
	scanf_s("%d", &Size);
	//sscanf_s(argv[1], "%d", &Size);
	//if (argc < 2) {
	//	printf("Error arguments command line");
	//	exit(1);
	//}

	double* MatrixA = NULL;
	if (!(MatrixA = (double*)malloc(sizeof(double)*Size*Size))) {
		printf("Error");
	}
	double* MatrixB = NULL;
	if (!(MatrixB = (double*)malloc(sizeof(double)*Size*Size))) {
		printf("Error");
	}
	double* MatrixC = NULL;
	if (!(MatrixC = (double*)malloc(sizeof(double)*Size*Size))) {
		printf("Error");
	}

	generator(MatrixA, MatrixB, MatrixC, Size);
	double start = omp_get_wtime();

	multiplication_block(MatrixA, MatrixB, MatrixC, Size);
	double end = omp_get_wtime();
	double delta = end - start;
	if (Size < 1000)
		print(MatrixA, MatrixB, MatrixC, Size);
	printf("Time for block parallel version: %.6f \n", delta);
	int s = 0;
	free(MatrixA);
	free(MatrixB);
	free(MatrixC);
	return s;
}