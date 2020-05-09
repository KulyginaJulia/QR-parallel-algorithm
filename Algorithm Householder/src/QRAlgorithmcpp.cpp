#include "../include/QRAlgorithm.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace QR;

QRAlgorithm::QRAlgorithm(int _n) {
	n = _n;
	A = new double[n * n];
	Q = new double[n * n];
	R = new double[n * n];
	P = new double[n * n]; // matrix of the v vectors
}

void QRAlgorithm::generator() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i*n + j] = rand() % 100;
		}
	}
}

double QRAlgorithm::Norma(double* x, int size) {
	double sum = 0;
	for (int i = 0; i < size; i++)
		sum += x[i] * x[i];
	return sum;
}

void QR::multiplication(double* MatrixA, double* MatrixB, double* &MatrixC, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			MatrixC[i*n + j] = 0;
			for (int k = 0; k < n; k++) {
				MatrixC[i * n + j] += MatrixA[i * n + k] * MatrixB[k * n + j];
			}
		}
	}
}

void QRAlgorithm::copy_matrix(double* &destiny, double* &source) {
	for (int i = 0; i < n*n; i++) {
		destiny[i] = source[i];
	}
}

double* QRAlgorithm::CreateMatrix_V(double* v) {
	double* V = new double[n * n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double val = v[i] * v[j];
			V[i* n + j] = val;
		}
	}

	return V;
}

void QRAlgorithm::check_result() {
	cout << "Check result matrices... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -7); //10 ^ (-7);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j) {
				if (abs(R[i * n + j]) > eps) {
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
	double* tmp_matrix = new double[n * n];

	multiplication(Q, R, tmp_matrix, n);
	// checking norma(A-QR)?
	double** MatrixA_QR = new double*[n];
	for (int i = 0; i < n; i++) {
		MatrixA_QR[i] = new double[n];
		for (int j = 0; j < n; j++) {
			MatrixA_QR[i][j] = A[i * n + j] - tmp_matrix[i * n + j];
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
	copy_matrix(R, A);	
	double* v = new double[n];
	double* x = new double[n]; //столбец матрицы
	for (int k = 0; k < n - 1; k++) {
		// cout << "Step " << k << endl;
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = R[j * n + k];
		}
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = R[k * n + k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i * n + k] + (i == k ? signum * norma_x : 0); 
			}

		}
		double norm = Norma(v, n);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i * n + k] = v[i] / sqrt(norm);
		}

		double* V = CreateMatrix_V(v);
		double* P_current = new double[n * n];
		double* tmpMatrix = new double[n * n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double val = (i == j ? 1 : 0) - V[i * n + j] * 2 / norm;
				P_current[i * n + j] = val;
			}
		}

		multiplication(P_current, R, tmpMatrix, n);
		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		copy_matrix(R, tmpMatrix);
		delete V;
		delete P_current;
		delete tmpMatrix;
	}
	delete v;
	delete x;
}

void PrimitiveQR::QSelector() {	
	double* v = new double[n];

	for (int k = n - 2; k >= 0; k--) {
		for (int i = 0; i < n; i++) {
			v[i] = P[i * n + k];
		}

		double* V = CreateMatrix_V(v);
		double* P_current = new double[n * n];
		double* tmpMatrix = new double[n * n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double val = (i == j ? 1 : 0) - V[i* n +j] * 2;
				P_current[i* n + j] = val;
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
		delete P_current;
		delete V;
	}
	delete v;
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

void RowHouseQR::row_house(double*& _A, double* v, int size, int k) {
	double norm = Norma(v, size);
	double beta = -2 / norm;

	double* w = new double[size];
	
	for (int i = k; i < size; i++) {
		w[i] = 0;
		for (int j = k; j < size; j++) {
			w[i] += _A[j* size + i] * v[j];
		}
		w[i] *= beta;
	}

	for (int i = k; i < size; i++) {
		for (int j = k; j < size; j++) {
			_A[i* size + j] = _A[i* size + j] + v[i] * w[j];
		}
	}

	delete[] w;
}


void QR::print_matrix(double* MatrixA, int Size) {
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			cout << fixed;
			cout.precision(6);
			cout << MatrixA[i* Size + j] << ' ';
		}
		cout << endl;
	}
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
				x[j] = R[j* n + k];
		}
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = R[k*n + k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i* n + k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v, n);

		// variant from GolubVanLoun
		row_house(R, v, n, 0);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i* n + k] = v[i] / sqrt(norm);
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

void MultiplicationQR::QSelector() {
	double* E = new double[n * n];

	// initialize E
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i* n + j] = 1;
			else
				E[i* n + j] = 0;
		}
	}
	copy_matrix(Q, E);

	for (int k = n - 2; k >= 0; k--) {
		double* v = new double[n];
		//cout << "Step = " << k << endl << endl;
		for (int i = 0; i < n; i++) {
			v[i] = P[i* n + k];
		}
		int r = n - k;

		double* Q_small = new double[r * r];
		for (int i = r - 1; i >= 0; i--) {
			for (int j = r - 1; j >= 0; j--)
				Q_small[i* (r) + j] = Q[(i + k) * (n) + (j + k)];
		}
		double* v_small = new double[r];

		for (int i = 0; i < r; i++) {
			v_small[i] = v[i + k];
		}

		row_house(Q_small, v_small, r, 0);

		for (int i = r - 1; i >= 0; i--) {
			for (int j = r - 1; j >= 0; j--)
				Q[(i + k) * (n) + (j + k)] = Q_small[i* (r) + j];
		}

		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
		//clear
		delete v;
		delete v_small;
		delete Q_small;
	}
}
void FinalSequenceQR::QRDecomposition() {
	copy_matrix(R, A);
	double* v = new double[n];
	double* x = new double[n]; //столбец матрицы
	for (int k = 0; k < n - 1; k++) {
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = R[j* n + k];
		}
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = R[k* n + k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i* n + k] + (i == k ? signum * norma_x : 0);
			}
		}
		double norm = Norma(v, n);

		// variant from GolubVanLoun
		row_house(R, v, n, 0);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i* n + k] = v[i] / sqrt(norm);
		}

		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;
	}
	delete v;
	delete x;
}

void FinalSequenceQR::QSelector() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q[i* n + j] = i == j ? 1 : 0;
		}
	}
	double* v = new double[n];
	for (int k = n - 2; k >= 0; k--) {
		//cout << "Step = " << k << endl << endl;
		for (int i = 0; i < n; i++) {
			v[i] = P[i* n + k];	
		}
	    row_house(Q, v, n, k);

		//if (n >= 100 && k % 10 == 0)
		//	cout << "... k = " << k << endl;

	}
	//clear
	delete v;
}

// Parallel
void ParallelQR::copy_matrix(double* &destiny, double* &source) {
#pragma omp parallel for
	for (int i = 0; i < n*n; i++) {
		destiny[i] = source[i];
	}
}

double ParallelQR::Norma(double* x, int size) {
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < size; i++)
		sum += x[i] * x[i];
	return sum;
}

void ParallelQR::row_house(double*& _A, double* v, int size, int k) {
	double norm = QRAlgorithm::Norma(v, size);
	double beta = -2 / norm;

	double* w = new double[size];
	int j = k;
#pragma omp parallel for
	for (int i = k; i < size; i++) {
		double ThreadResult = 0;
#pragma omp parallel for reduction(+:ThreadResult)
		for (int j = k; j < size; j++)
			ThreadResult += _A[j*size + i] * v[j];
		w[i] = ThreadResult;
	}

//#pragma omp parallel for private(j)
//	for (int i = k; i < size; i++) {
//		w[i] = 0;
//		for (j = k; j < size; j++) {
//			w[i] += _A[j* size + i] * v[j];
//		}
//		//w[i] *= beta;
//	}

#pragma omp parallel for private(j)
	for (int i = k; i < size; i++) {
		for (j = k; j < size; j++) {
			_A[i* size + j] += v[i] * w[j] * beta;
		}
	}

	delete[] w;
}

void ParallelQR::QRDecomposition() {
	copy_matrix(R, A);
	double* v = new double[n];
	double* x = new double[n]; //столбец матрицы

	for (int k = 0; k < n - 1; k++) {

//#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			if (j < k)
				x[j] = 0;
			else
				x[j] = R[j* n + k];
		}
		double norma_x = sqrt(Norma(x, n));

//#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			int signum = R[k* n + k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = R[i* n + k] + (i == k ? signum * norma_x : 0);
			}
		}
		double norm = Norma(v, n);

		// variant from GolubVanLoun
		row_house(R, v, n, 0);

//#pragma omp parallel for
		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i* n + k] = v[i] / sqrt(norm);
		}
	}
	delete v;
	delete x;
}

void ParallelQR::QSelector() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q[i* n + j] = i == j ? 1 : 0;
		}
	}
	double* v = new double[n];
	for (int k = n - 2; k >= 0; k--) {
//#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			v[i] = P[i* n + k];
		}

		row_house(Q, v, n, k);

	}
	//clear
	delete v;
}

void GivensRotation::givens(double a, double b, double &c, double &s) {
	if (b == 0) {
		c = 1;
		s = 0;
		return;
	}
	//double r = abs(sqrt(a*a + b*b));
	//c = a / r;
	//s = -b / r;
	if (abs(b) > abs(a)) {
		double tao = -a / b;
		s = 1 / sqrt(1 + tao * tao); // -1?
		c = s * tao;
	}
	else {
        //double signum_a = a > 0 ? 1 : -1;
		double tao = -b / a;
		c = 1 / sqrt(1 + tao * tao);
		s = c * tao;
	}
}
// A = GT*A
void GivensRotation::rowRotation(double*& A, double c, double s, int _i, int _j) {

	for (int j = _j; j < n; j++) {
		int index1 = (_i - 1) * n + j;
		int index2 = _i * n + j;
		double tao1 = A[index1];
		double tao2 = A[index2];
		A[index1] = c * tao1 - s * tao2;
		A[index2] = s * tao1 + c * tao2;
	}
}
// A = A*GT
void GivensRotation::columnRotation(double*& A, double c, double s, int _i, int _j) {

	for (int i = 0; i < n; i++) {
		int index1 = i * n + _j;
		int index2 = i * n + _j + 1;
		double tao1 = A[index1];
		double tao2 = A[index2];
		A[index1] = c * tao1 - s * tao2;
		A[index2] = s * tao1 + c * tao2;
	}
}

int GivensRotation::calculateP(double c, double s) {
	if (c == 0) {
		return 1;
	}
	else if (abs(s) < abs(c)) {
		int signum = c > 0 ? 1 : -1;
		return signum * s / 2;
	}
	else {
		int signum = s > 0 ? 1 : -1;
		return 2 * signum / c;
	}
}

void GivensRotation::QRDecomposition() {
	copy_matrix(R, A);

	cout << "Matrix R = " << endl;
	print_matrix(R, n);
	int count = 0;
	for (int j = 0; j < n - 1; j++) {
		for (int i = n-1; i > 0; i -= j + 1) {
			cout << "Step i = " << i << " j = " << j <<  endl;
			double c, s;
			givens(R[(i - 1) * n + j], R[i * n + j], c, s);
			cout << "c = " << c << "  s = " << s << endl;
			rowRotation(R, c, s, i, j);
			//p[i * (n) + j] = calculateP(c, s);
			//cout << "index = " << i * (n) + j << "  p = " << p[i * (n) + j] << endl;
			cout << "index = " << count << endl;
			coeffsC[count] = c;
			coeffsS[count] = s;
			//print_matrix(R, n);
			count++;
		}
	}

	//print_matrix(R, n);
}
void GivensRotation::calculateCSfromP(double p, double &c, double &s) {
	if (p == 1) {
		c = 0;
		s = 1;
	}
	else if (abs(p) < 1) {
		s = 2 * p;
		c = sqrt(1 - s * s);
	}
	else {
		c = 2 / p;
		s = sqrt(1 - c * c);
	}
}

void copyG1toQ(double*& Q, double c, double s, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int index = i * n + j;
			Q[index] = i == j ? 1 : 0;
		}
	}
	int index = (n - 2) * n + n - 2;
	Q[(n - 2) * n + n - 2] = c;
	Q[(n - 2) * n + n - 1] = -s;
	Q[(n - 1) * n + n - 2] = s;
	Q[(n - 1) * n + n - 1] = c;

}
void GivensRotation::QSelector() {
	int count = 0;
	for (int j = 0; j < n - 1; j++) {
		for (int i = n - 1; i > 0; i -= j + 1) {
			double c = coeffsC[count], s = coeffsS[count];
			count++;
			//calculateCSfromP(p[i* (n) + j], c, s);
			if (count ==  1) {
				copyG1toQ(Q, c, s, n);
				//cout << "G1: " << endl;
				//print_matrix(Q, n);
				continue;
			}
			columnRotation(Q, c, s, i, j);
			//cout << "After rotation " << endl;
			//print_matrix(Q, n);
		}
	}

	//print_matrix(Q, n);
}
