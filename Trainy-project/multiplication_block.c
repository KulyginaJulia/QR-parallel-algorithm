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
void multiplication_block(double* MatrixA, double* MatrixB, double* MatrixC, int N, int ThreadNum) {
	double sqThreadNum = sqrt((double)ThreadNum);
	int GridSize = (int)(sqThreadNum);
	int BlockSize = N / GridSize;
	omp_set_num_threads(ThreadNum);
#pragma omp parallel
	{ 		
		int ThreadID = omp_get_thread_num();
		int RowIndex = ThreadID / GridSize;
		int ColIndex = ThreadID % GridSize;
		for (int iter = 0; iter < GridSize; iter++) {
			for (int i = RowIndex * BlockSize; i < (RowIndex + 1) * BlockSize; i++)
				for (int j = ColIndex * BlockSize; j < (ColIndex + 1) * BlockSize; j++)
					for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++)
						MatrixC[i*N + j] += MatrixA[i*N + k] * MatrixB[k*N + j];
		} 

	}// pragma omp parallel
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
	int ThredNum = 0;
	scanf_s("%d", &Size);
	scanf_s("%d", &ThredNum);
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

	multiplication_block(MatrixA, MatrixB, MatrixC, Size, ThredNum);
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