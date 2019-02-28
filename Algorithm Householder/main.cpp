#include <iostream>
#include <list>
#include <vector>
#include <cmath>

using namespace std;

double Normirovka(double *s, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += abs(s[i]) * abs(s[i]);
	return sum;
}

double** CreateMatrix_V(double* v, int size) {
	double** V = new double*[size];
	for (int i = 0; i < size; i++) {
		V[i] = new double[size];
		for (int j = 0; j < size; j++) {
			double val = v[i] * v[j];
			V[i][j] = val;
		}
	}

	return V;
}

void print_matrix(double** &MatrixA, int Size) {
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << MatrixA[i][j] << ' ';
		}
		cout << endl;
	}
}

void copy_matrix(double** &destiny, double** &source, int Size) {
	destiny = new double*[Size];
	for (int i = 0; i < Size; i++) {
		destiny[i] = new double[Size];
		for (int j = 0; j < Size; j++) {
			destiny[i][j] = source[i][j];
		}
	}
}

double** multiplication(double** MatrixA, double** MatrixB, int Size) {
	double** MatrixC = new double*[Size];
	for (int i = 0; i < Size; i++)
		MatrixC[i] = new double[Size];

	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			MatrixC[i][j] = 0;
			for (int k = 0; k < Size; k++) {
				MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
			}
		}
	}
	//for (int i = 0; i < Size; i++) {
	//	for (int j = 0; j < Size; j++) {
	//		//MatrixC[j].push_back(0);
	//		//floor(atof(buf) * 100 + 0.5) / 100
	//		MatrixC[i][j] = std::floor(MatrixC[i][j] * 100 + 0.5) / 100;
	//	}
	//}
	return MatrixC;
}

void check_result(double** MatrixA, double** MatrixQ, double** MatrixR, int size) {
	cout << "Check result matrices... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i > j)
				if (MatrixR[i][j] != 0) {
					flag_ok = false;
					cout << "R[" << i << "][" << j << "] = " << MatrixR[i][j] << "!= 0" << endl;
				}
		}
	}
	if (flag_ok)
		cout << "  OK" << endl;
	else
		cout << "  ERR" << endl;
	cout << "Done." << endl;
	flag_ok = true;

	cout << "Check A = Q*R: " << endl;
	double** tmp_matrix = multiplication(MatrixQ, MatrixR, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (MatrixA[i][j] != round(tmp_matrix[i][j])) {
				flag_ok = false;
				cout << "Mismatch! A[" << i << "][" << j << "] = " << MatrixA[i][j] << "!= " << round(tmp_matrix[i][j]) << " = Q*R[" << i << "][" << j << "]" << endl;

			}
		}
	}

	if (flag_ok)
		cout << "  OK" << endl;
	else
		cout << "  ERR" << endl;
	cout << "Done." << endl;
}

int main() {
	double E[3][3] = { {1, 0, 0},
				   {0, 1, 0},
				   {0, 0, 1 }
	};
	int n = 3;
	double** A = new double*[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
	}	
	A[0][0] = 1;	A[0][1] = -2;	A[0][2] = 1;
	A[1][0] = 2; 	A[1][1] = -1; 	A[1][2] = 3;
	A[2][0] = 2; 	A[2][1] = 3; 	A[2][2] = -4;
	double** A_ = A;
	double** Q;

	for (int k = 0; k < n - 1; k++) {
		cout << "Step " << k << endl;
		double* v = new double[n];
		double* x = new double[n]; //столбец матрицы
		cout << "  x = " << endl;
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else 
				x[j] = A_[j][k];
			cout << x[j] << endl;
		}
		cout << endl;
		double norma_x = sqrt(Normirovka(x, n));
		cout << "  Norma x = " << norma_x << endl;
		for (int i = 0; i < n; i++) {
			int signum = A_[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = A_[i][k] + signum * norma_x * E[i][k];
			}

		}
		cout << "  v = ";
		for (int j = 0; j < n; j++) {
			cout << v[j] << " ";
		}

		double norm = Normirovka(v, n);
		cout << "  Quad Norma v = " << norm << endl;

		// double sqrt_norm = sqrt(norm);
		double** V = CreateMatrix_V(v, n);
		cout << "  Matrix V:" << endl;
		print_matrix(V, n);

		double** P_current = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2 / norm;
				P_current[i][j] = val;
			}
		}
		cout << "  Matrix P:" << endl;
		print_matrix(P_current, n);

		
		A_ = multiplication(P_current, A_, n);

		cout << "  Matrix R:" << endl;
		print_matrix(A_, n);
		int g = 0;
		if (k == 0)
			copy_matrix(Q, P_current, n);
		else
			Q = multiplication(P_current, Q, n);
		// clear
		for (int i = 0; i < n; i++)
			delete P_current[i];

		delete P_current;
		delete v;
	}
	double** R = A_;

	check_result(A, Q, R, n);
	return 0;
}