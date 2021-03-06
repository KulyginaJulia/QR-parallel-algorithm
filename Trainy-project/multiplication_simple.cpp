#include <omp.h>
#include <iostream>

using namespace std;

void generator(double* &MatrixA, double* &MatrixB, double* &MatrixC, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrixA[i*N + j] = rand() % 100;
			MatrixB[i*N + j] = rand() % 100;
			MatrixC[i*N + j] = 0;
		}
	}
}
void multiplication(double* &MatrixA, double* &MatrixB, double* &MatrixC, int N) {
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<N; j++) {
			MatrixC[i*N + j] = 0;
			for (int k = 0; k<N; k++) {
				MatrixC[i*N + j] = MatrixC[i*N + j] + MatrixA[i*N + k] * MatrixB[k*N + j];
			}
		}
	}
}
void print(double* &MatrixA, double* &MatrixB, double* &MatrixC, int Size) {
	cout << "Matrix A:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixA[i*Size + j] << ' ';
		}
		cout << endl;
	}

	cout << "Matrix B:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixB[i*Size + j] << ' ';
		}
		cout << endl;
	}
	cout << "Matrix C:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixC[i*Size + j] << ' ';
		}
		cout << endl;
	}
}
int main() {
	
	int Size;
	std::cin >> Size;

	double* MatrixA = new double[Size*Size];
	double* MatrixB = new double[Size*Size];
	double* MatrixC = new double[Size*Size];


	generator(MatrixA, MatrixB, MatrixC, Size);
	double start = omp_get_wtime();
	multiplication(MatrixA, MatrixB, MatrixC, Size);
	double end = omp_get_wtime();
	double delta = end - start;
	if (Size < 1000)
		print(MatrixA, MatrixB, MatrixC, Size);
	cout << "Time for simple version: " << delta << endl;
	int s = 0;
}