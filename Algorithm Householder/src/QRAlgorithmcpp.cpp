#include "../include/QRAlgorithm.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace QR;

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

double QRAlgorithm::Norma(double* x, int size) {
	double sum = 0;
	for (int i = 0; i < size; i++)
		sum += x[i] * x[i];
	return sum;
}

void QR::multiplication(double** MatrixA, double** MatrixB, double** &MatrixC, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			MatrixC[i][j] = 0;
			for (int k = 0; k < n; k++) {
				MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
			}
		}
	}
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
			if (i > j) {
				if (abs(R[i][j]) > eps) {
					flag_ok = false;
					//cout << "R[" << i << "][" << j << "] = " << MatrixR[i][j] << "!= 0" << endl;
					break;
				}
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
	double** tmp_matrix = new double*[n];
	for (int i = 0; i < n; i++) {
		tmp_matrix[i] = new double[n];
	}

	multiplication(Q, R, tmp_matrix, n);
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
	//double** A_ = new double*[n];

	// initialize E
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
		//A_[i] = new double[n];
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}
	copy_matrix(R, A);
	for (int k = 0; k < n - 1; k++) {
		// cout << "Step " << k << endl;
		double* v = new double[n];
		double* x = new double[n]; //столбец матрицы
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = R[j][k];
		}
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = R[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i][k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v, n);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i][k] = v[i] / sqrt(norm);
		}

		double** V = CreateMatrix_V(v);
		double** P_current = new double*[n];
		double** tmpMatrix = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			tmpMatrix[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2 / norm;
				P_current[i][j] = val;
			}
		}

		multiplication(P_current, R, tmpMatrix, n);
		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		copy_matrix(R, tmpMatrix);
		for (int i = 0; i < n; i++) {
			delete P_current[i];
			delete V[i];
			delete tmpMatrix[i];
		}
		delete V;
		delete P_current;
		delete v;
		delete tmpMatrix;
	}
	//R = A_;
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
		double** tmpMatrix = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			tmpMatrix[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2;
				P_current[i][j] = val;
			}
		}

		if (k == n - 2) {
			copy_matrix(Q, P_current);
		}
		else {
			multiplication(P_current, Q, tmpMatrix, n);
			copy_matrix(Q, tmpMatrix);
		}
		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		// clear
		for (int i = 0; i < n; i++) {
			delete V[i];
			delete P_current[i];
			delete tmpMatrix[i];
		}

		delete P_current;
		delete V;
		delete v;
	}
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

void RowHouseQR::row_house(double**& _A, double* v, int size) {
	double norm = Norma(v, size);
	double beta = -2 / norm;

	double* w = new double[size];
	
	for (int i = 0; i < size; i++) {
		w[i] = 0;
		for (int j = 0; j < size; j++) {
			w[i] += _A[j][i] * v[j];
		}
		w[i] *= beta;
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			_A[i][j] = _A[i][j] + v[i] * w[j];
		}
	}

	delete[] w;
}

void RowHouseQR::QRDecomposition() {
	double** E = new double*[n];

	// initialize E
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}

	copy_matrix(R, A);

	for (int k = 0; k < n - 1; k++) {
		double* v = new double[n];
		double* x = new double[n]; //столбец матрицы
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = R[j][k];
		}
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = R[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i][k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v, n);

		// variant from GolubVanLoun
		row_house(R, v, n);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i][k] = v[i] / sqrt(norm);
		}

		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		delete v;
	}

	for (int i = 0; i < n; i++) {
		delete E[i];
	}
	delete E;
}

void QR::print_matrix(double** MatrixA, int Size) {
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			//cout << fixed;
			//cout.precision(6);
			cout << MatrixA[i][j] << ' ';
		}
		cout << endl;
	}
}


void MultiplicationQR::QSelector() {
	double** E = new double*[n];

	// initialize E
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}
	copy_matrix(Q, E);

	for (int k = n - 2; k >= 0; k--) {
		double* v = new double[n];
		//cout << "Step = " << k << endl << endl;
		for (int i = 0; i < n; i++) {
			v[i] = P[i][k];	
		}
		int r = n - k;

		double** Q_small = new double*[r];
		for (int i = r-1; i >= 0; i--) {
			Q_small[i] = new double[r];
			for (int j = r-1; j >= 0; j--)
				Q_small[i][j] = Q[i+k][j+k];
		}
		double* v_small = new double[r];

		for (int i = 0; i < r; i++) {
			v_small[i] = v[i + k];
		}

	    row_house(Q_small, v_small, r);

		for (int i = r-1; i >= 0; i--) {
			for (int j = r-1; j >= 0; j--)
				Q[i + k][j + k] = Q_small[i][j];
		}

		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		//clear
		delete v;
		delete v_small;
		for (int i = 0; i < r; i++)
			delete[] Q_small[i];
		delete Q_small;
	}
}