#include <omp.h>
#include <iostream>

using namespace std;

void generator(double** MatrixA, double** MatrixB, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrixA[i][j] = rand() % 10;
			MatrixB[i][j] = rand() % 10;
		}
	}
}
void multiplication(double** MatrixA, double** MatrixB, double** MatrixC, int N) {
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<N; j++) {
			MatrixC[i][j] = 0;
			for (int k = 0; k<N; k++) {
				MatrixC[i][j] = MatrixC[i][j] + MatrixA[i][k] * MatrixB[k][j];
			}
		}
	}
}
void print(double** MatrixA, double** MatrixB, double** MatrixC, int Size) {
	cout << "Matrix A:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixA[i][j] << ' ';
		}
		cout << endl;
	}

	cout << "Matrix B:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixB[i][j] << ' ';
		}
		cout << endl;
	}
	cout << "Matrix C:" << endl;
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixC[i][j] << ' ';
		}
		cout << endl;
	}
}
int main() {
	
	int Size;
	std::cin >> Size;

	double** MatrixA = new double*[Size];
	double** MatrixB = new double*[Size];
	double** MatrixC = new double*[Size];

	for (int i = 0; i < Size; i++) {
		MatrixA[i] = new double[Size];
		MatrixB[i] = new double[Size];
		MatrixC[i] = new double[Size];
	}

	generator(MatrixA, MatrixB, Size);
	double start = omp_get_wtime();
	multiplication(MatrixA, MatrixB, MatrixC, Size);
	double end = omp_get_wtime();
	double delta = end - start;
	print(MatrixA, MatrixB, MatrixC, Size);
	cout << "Time for simple version: " << delta << endl;
	int s = 0;
}