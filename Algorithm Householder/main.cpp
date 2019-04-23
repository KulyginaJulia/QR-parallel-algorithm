#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>

using namespace std;
//using namespace Eigen;


void generator(double** &MatrixA, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrixA[i][j] = rand() % 100;
		}
	}
}

double Norma(double *s, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += abs(s[i]) * abs(s[i]);
	return sum;
}

void Normirovka(double* &s, int n, double norma) {
	for (int i = 0; i < n; i++)
		s[i] = s[i] / sqrt(norma);
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
			//cout << fixed;
			//cout.precision(6);
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
	return MatrixC;
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

void check_result(double** &MatrixA, double** &MatrixQ, double** &MatrixR, int size) {
	cout << "Check result matrices... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -14); //10 ^ (-14);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i > j)
				if (abs(MatrixR[i][j]) > eps) {
					flag_ok = false;
					//cout << "R[" << i << "][" << j << "] = " << MatrixR[i][j] << "!= 0" << endl;
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
	// checking norma(A-QR)?
	double** MatrixA_QR = new double*[size];
	for (int i = 0; i < size; i++) {
		MatrixA_QR[i] = new double[size];
		for (int j = 0; j < size; j++) {
			MatrixA_QR[i][j] = MatrixA[i][j] - tmp_matrix[i][j];
		}
	}
	//cout << "Matrix A-QR" << endl;
	//print_matrix(MatrixA_QR, size);
	double *vector_column = new double[size];
	
	for (int i = 0; i < size; i++)
	{
		vector_column[i] = 0;
		for (int j = 0; j < size; j++)
			vector_column[i] += abs(MatrixA_QR[i][j]);
	}
	double max_abs = vector_column[0];
	for (int i = 0; i < size; i++)
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

void check_with_eigen_result(double** &MatrixQ, double** &MatrixR, Eigen::MatrixXd MatrixQE, Eigen::MatrixXd MatrixRE, int size) {
	cout << "Check result matrices with eigen results... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -14); //10 ^ (-14);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i < j)
				if (abs(MatrixR[i][j] - MatrixRE(i, j)) > eps) {
					flag_ok = false;
					//cout << "R[" << i << "][" << j << "] = " << MatrixR[i][j] << "!= R`["
					//	<< i << "][" << j << "] = " << MatrixRE(i, j) << endl;
				}
		}
	}
	if (flag_ok)
		cout << "  OK" << endl;
	else
		cout << "  ERR" << endl;
	cout << "Done." << endl;
	flag_ok = true;

	cout << "Check Q: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (abs(MatrixQ[i][j] - MatrixQE(i, j)) > eps) {
				flag_ok = false;
				//cout << "Q[" << i << "][" << j << "] = " << MatrixQ[i][j] << "!= Q`["
				//	<< i << "][" << j << "] = " << MatrixQE(i, j) << endl;
			}
		}
	}
	if (flag_ok)
		cout << "  OK" << endl;
	else
		cout << "  ERR" << endl;
	cout << "Done." << endl;
}

double** calculate_Q_matrix(double** &P, int n, double** &E, double** &A_) {

	double** Q = new double*[n];
	for (int i = 0; i < n; i++)
		Q[i] = new double[n];

	for (int k = n - 2; k >= 0; k--) {
		double* v = new double[n];
		for (int i = 0; i < n; i++) {
			v[i] = P[i][k];
		}

		double** V = CreateMatrix_V(v, n);
		double** P_current = new double*[n];
		for (int i = 0; i < n; i++) {
			P_current[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - V[i][j] * 2;
				P_current[i][j] = val;
			}
		}

		if (k == n - 2)
			copy_matrix(Q, P_current, n);
		else
			Q = multiplication(P_current, Q, n);
		// clear
		for (int i = 0; i < n; i++)
			delete P_current[i];

		delete P_current;
		delete v;
	}
	return Q;
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

double** row_house(double** A, double* v, int n) {
	double norm = Norma(v, n);
	double beta = -2 / norm;
	double** At = Transpose_(A, n);
	double* w = multiplication_matrix_on_vector(At, v, n);
	for (int i = 0; i < n; i++)
		w[i] *= beta;
	double** v_wt = multiplication_column_on_row(v, w, n);

	double** A_ = new double*[n];
	for (int i = 0; i < n; i++) {
		A_[i] = new double[n];
		for (int j = 0; j < n; j++) {
			A_[i][j] = A[i][j] + v_wt[i][j];
		}
	}
	return A_;
}


void QRDecomposition(int n, double** A_, double** &Q, double** &R) {
	double** E = new double*[n];
	for (int i = 0; i < n; i++) {
		E[i] = new double[n];
	}
	// initialize E
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
	double** P = new double*[n]; // matrix of the v vectors
	for (int i = 0; i < n; i++) {
		P[i] = new double[n];
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
		double norma_x = sqrt(Norma(x, n));
		for (int i = 0; i < n; i++) {
			int signum = A_[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v[i] = 0;
			}
			else {
				v[i] = A_[i][k] + signum * norma_x * E[i][k];
			}

		}
		double norm = Norma(v, n);
		//Normirovka(v, n, norm);

		// variant from GolubVanLoun
		A_ = row_house(A_, v, n);

		// add v to P matrix as a column
		for (int i = 0; i < n; i++) {
			P[i][k] = v[i] / sqrt(norm);
		}

		//double** V = CreateMatrix_V(v, n);
		//double** P_current = new double*[n];
		//for (int i = 0; i < n; i++) {
		//	P_current[i] = new double[n];
		//	for (int j = 0; j < n; j++) {
		//		double val = E[i][j] - V[i][j] * 2 / norm;
		//		P_current[i][j] = val;
		//	}
		//}
		// in lob
		//A_ = multiplication(P_current, A_, n);

		delete v;
	}
	Q = calculate_Q_matrix(P, n, E, A_);
	R = A_;
}



int main(int argc, char **argv) {
	int n = atoi(argv[1]);
	int mode = atoi(argv[2]);
	double** A = new double*[n];

	double** Q, **R;
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
	}
	cout.scientific;
	if (mode == 0) {
		n = 3;
		A[0][0] = 1;	A[0][1] = -2;	A[0][2] = 1;
		A[1][0] = 2; 	A[1][1] = -1; 	A[1][2] = 3;
		A[2][0] = 2; 	A[2][1] = 3; 	A[2][2] = -4;
	}
	else 
		generator(A, n);

	double** A_ = A;
	double start = omp_get_wtime();


	QRDecomposition(n, A_, Q, R);

	double end = omp_get_wtime();
	double delta = end - start;
	// Q = Transpose_(Q, n);
	//std::cout << "Matrix Q:" << std::endl;
	//print_matrix(Q, n);
	//std::cout << "Matrix R:" << std::endl;
	//print_matrix(R, n);

	check_result(A, Q, R, n);
	std::cout << "Time for simple version: " << delta << std::endl;


	Eigen::MatrixXd m(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			m(i, j) = A[i][j];
	}
	auto qr = m.householderQr(); // casted type Matrix -> Householder matrix

	start = omp_get_wtime();
	
	auto r = qr.matrixQR(); // qr decomposition, r - upper triangle, vectors v - lower triangle

	end = omp_get_wtime();

	auto q = qr.householderQ(); // selection q
	Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(n, n));
	thinQ = q * thinQ;

	delta = end - start;
	cout << "Time for eigen version: " << delta << endl;
	//cout << "Matrix q " << endl << thinQ << endl;
	//cout << "Matrix r " << endl;
	//cout << r << endl;

	check_with_eigen_result(Q, R, thinQ, r, n);

	return 0;
}