#include "../include/QRAlgorithm.h"
#include <iostream>
#include <cmath>

using namespace std;

QRAlgorithm::QRAlgorithm(int _n) {
	n = _n;
	A = new double*[n];
	Q = new double*[n];
	R = new double*[n];
	P = new double*[n]; // matrix of the v vectors
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
		Q[i] = new double[n];
		R[i] = new double[n];
		P[i] = new double[n];
	}
}

void QRAlgorithm::generator() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = rand() % 100;
		}
	}
}

double QRAlgorithm::Norma(double* x) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += abs(x[i]) * abs(x[i]);
	return sum;
}

double** QRAlgorithm::multiplication(double** MatrixA, double** MatrixB) {
	double** MatrixC = new double*[n];
	for (int i = 0; i < n; i++)
		MatrixC[i] = new double[n];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			MatrixC[i][j] = 0;
			for (int k = 0; k < n; k++) {
				MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
			}
		}
	}
	return MatrixC;
}

void QRAlgorithm::copy_matrix(double** &destiny, double** &source) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			destiny[i][j] = source[i][j];
		}
	}
}

double** QRAlgorithm::CreateMatrix_V(double* v) {
	double** V = new double*[n];
	for (int i = 0; i < n; i++) {
		V[i] = new double[n];
		for (int j = 0; j < n; j++) {
			double val = v[i] * v[j];
			V[i][j] = val;
		}
	}

	return V;
}

void QRAlgorithm::check_result() {
	cout << "Check result matrices... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -10); //10 ^ (-10);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j)
				if (abs(R[i][j]) > eps) {
					flag_ok = false;
					//cout << "R[" << i << "][" << j << "] = " << MatrixR[i][j] << "!= 0" << endl;
					break;
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
	double** tmp_matrix = multiplication(Q, R);
	// checking norma(A-QR)?
	double** MatrixA_QR = new double*[n];
	for (int i = 0; i < n; i++) {
		MatrixA_QR[i] = new double[n];
		for (int j = 0; j < n; j++) {
			MatrixA_QR[i][j] = A[i][j] - tmp_matrix[i][j];
		}
	}
	//cout << "Matrix A-QR" << endl;
	//print_matrix(MatrixA_QR, size);
	double *vector_column = new double[n];

	for (int i = 0; i < n; i++)
	{
		vector_column[i] = 0;
		for (int j = 0; j < n; j++)
			vector_column[i] += abs(MatrixA_QR[i][j]);
	}
	double max_abs = vector_column[0];
	for (int i = 0; i < n; i++)
	{
		//cout << "vector_column = " << vector_column[i] << endl;
		if (vector_column[i] > max_abs)
			max_abs = vector_column[i];
	}
	cout << "max abs = " << max_abs << endl;
	if (max_abs < eps)
		cout << "  OK" << endl;
	else
		cout << "  ERR" << endl;
	cout << "Done." << endl;
}

void PrimitiveQR::QRDecomposition() {
	double** E = new double*[n];
	double** A_ = new double*[n];
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
		A_[i] = new double[n];
	}
	// initialize E
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
			A_[i] = A[i];
		}

	for (int k = 0; k < n - 1; k++) {
		// cout << "Step " << k << endl;
		double* v = new double[n];
		double* x = new double[n]; //столбец матрицы
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = A_[j][k];
		}
		double norma_x = sqrt(Norma(x));
		for (int i = 0; i < n; i++) {
			int signum = A_[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = A_[i][k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i][k] = v[i] / sqrt(norm);
		}

		double** V = CreateMatrix_V(v);
		double** P_current = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2 / norm;
				P_current[i][j] = val;
			}
		}
		A_ = multiplication(P_current, A_);
		if (n >= 100 && k % 10 == 0)
			cout << "... k = " << k << endl;
		delete v;
	}
	//Q = calculate_Q_matrix(P, n, E, A_);
	R = A_;
}

void PrimitiveQR::QSelector() {
	double** E = new double*[n];
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
	}
	// initialize E
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}

	for (int k = n - 2; k >= 0; k--) {
		double* v = new double[n];
		for (int i = 0; i < n; i++) {
			v[i] = P[i][k];
		}

		double** V = CreateMatrix_V(v);
		double** P_current = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2;
				P_current[i][j] = val;
			}
		}

		if (k == n - 2)
			copy_matrix(Q, P_current);
		else
			Q = multiplication(P_current, Q);
		// clear
		for (int i = 0; i < n; i++)
			delete P_current[i];

		delete P_current;
		for (int i = 0; i < n; i++)
			delete V[i];

		delete V;
		delete v;
	}
}

double** Transpose_(double** &MatrixQ, int size) {
	double** MatrixQ_ = new double*[size];
	double tmp = 0;
	for (int i = 0; i < size; i++) {
		MatrixQ_[i] = new double[size];
		for (int j = 0; j < size; j++) {
			MatrixQ_[i][j] = MatrixQ[j][i];
		}
	}

	return MatrixQ_;
}

double* multiplication_matrix_on_vector(double** MatrixA, double* VectorB, int Size) {
	double* VectorC = new double[Size];

	for (int i = 0; i < Size; i++) {
		VectorC[i] = 0;
		for (int j = 0; j < Size; j++) {
			VectorC[i] += MatrixA[i][j] * VectorB[j];
		}
	}
	return VectorC;
}

double** multiplication_column_on_row(double* VectorA, double* VectorB, int Size) {
	double** MatrixC = new double*[Size];
	for (int i = 0; i < Size; i++)
		MatrixC[i] = new double[Size];

	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			MatrixC[i][j] = VectorA[i] * VectorB[j];
		}
	}
	return MatrixC;
}

double** RowHouseQR::row_house(double** _A, double* v) {
	double norm = Norma(v);
	double beta = -2 / norm;
	double** At = Transpose_(_A, n);
	double* w = multiplication_matrix_on_vector(At, v, n);
	for (int i = 0; i < n; i++)
		w[i] *= beta;
	double** v_wt = multiplication_column_on_row(v, w, n);

	double** A_ = new double*[n];
	for (int i = 0; i < n; i++) {
		A_[i] = new double[n];
		for (int j = 0; j < n; j++) {
			A_[i][j] = _A[i][j] + v_wt[i][j];
		}
	}
	for (int i = 0; i < n; i++)
		delete[] At[i];
	delete[] At;

	delete[] w;
	for (int i = 0; i < n; i++)
		delete[] v_wt[i];
	delete[] v_wt;
	return A_;
}

void RowHouseQR::QRDecomposition() {
	double** E = new double*[n];
	double** A_ = new double*[n];
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
		A_[i] = new double[n];
	}
	// initialize E
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
			A_[i] = A[i];
		}

	for (int k = 0; k < n - 1; k++) {
		// cout << "Step " << k << endl;
		double* v = new double[n];
		double* x = new double[n]; //столбец матрицы
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = A_[j][k];
		}
		double norma_x = sqrt(Norma(x));
		for (int i = 0; i < n; i++) {
			int signum = A_[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = A_[i][k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v);

		// variant from GolubVanLoun
		A_ = row_house(A_, v);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i][k] = v[i] / sqrt(norm);
		}

		if (n >= 100 && k % 10 == 0)
			cout << "... k = " << k << endl;
		delete v;
	}
	R = A_;

}