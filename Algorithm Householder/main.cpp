#include <iostream>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>
#include "include/QRAlgorithm.h"

using namespace std;
using namespace QR;
//using namespace Eigen;

void check_with_eigen_result(double** MatrixQ, double** MatrixR, Eigen::MatrixXd MatrixQE, Eigen::MatrixXd MatrixRE, int size) {
	cout << "Check result matrices with eigen results... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -10); //10 ^ (-10);
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
				break;
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
	cout << "Please choose mode: 0 - Test mode with definitly 3x3 matrix, 1 - primitive version, 2 - rowHouse version, 3 - multipilication version" << endl;

	int n, mode;
	cin >> mode;
	cout << "Please enter size of matrix: ";
	cin >> n;
	cout << endl;
	//mode = 3;
	//n = 1000;
	cout.scientific;
	QRAlgorithm *qrAlgo;
	PrimitiveQR *primitiveQR = new PrimitiveQR(n);
	MultiplicationQR *multiplicationQR = new MultiplicationQR(n);
	RowHouseQR *rowHouseQR = new RowHouseQR(n);

	switch (mode) {
	case 0:
	{
		//qrAlgo = new MultiplicationQR(n);
		qrAlgo = multiplicationQR;
		qrAlgo->A[0][0] = 1;	qrAlgo->A[0][1] = -2;	qrAlgo->A[0][2] = 1;
		qrAlgo->A[1][0] = 2; 	qrAlgo->A[1][1] = -1; 	qrAlgo->A[1][2] = 3;
		qrAlgo->A[2][0] = 2; 	qrAlgo->A[2][1] = 3; 	qrAlgo->A[2][2] = -4;
		break;
	}
	case 1: {
		//qrAlgo = new PrimitiveQR(n);
		qrAlgo = primitiveQR;
		qrAlgo->generator();
		break;
	}
	case 2:
	{
		//qrAlgo = new RowHouseQR(n);
		qrAlgo = rowHouseQR;
		qrAlgo->generator();
		break;
	}
	case 3: 
	{
		//qrAlgo = new MultiplicationQR(n);
		qrAlgo = multiplicationQR;
		qrAlgo->generator();
		break;
	}
	default:
		break;

	}

	double start_decomp = omp_get_wtime();

	qrAlgo->QRDecomposition();

	double end_decomp = omp_get_wtime();
	double delta_decomp = end_decomp - start_decomp;



	cout << "Done QR decomposition" << endl;
	std::cout << "Time for simple version of decomposition: " << delta_decomp << std::endl;
	double start_q = omp_get_wtime();

	qrAlgo->QSelector();

	double end_q = omp_get_wtime();
	double delta_q = end_q - start_q;
	std::cout << "Time for simple version of Q: " << delta_q << std::endl;

	//std::cout << "Matrix Q:" << std::endl;
	//print_matrix(qrAlgo->get_Q(), n);

	//std::cout << "Matrix R:" << std::endl;
	//print_matrix(qrAlgo->get_R(), n);

	qrAlgo->check_result();

	Eigen::MatrixXd m(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			m(i, j) = qrAlgo->A[i][j];
	}
	auto qr = m.householderQr(); // casted type Matrix -> Householder matrix

	start_decomp = omp_get_wtime();
	
	auto r = qr.matrixQR(); // qr decomposition, r - upper triangle, vectors v - lower triangle

	end_decomp = omp_get_wtime();
	delta_decomp = end_decomp - start_decomp;
	cout << "Time for eigen version of decomposition: " << delta_decomp << endl;

	start_q = omp_get_wtime();
	auto q = qr.householderQ(); // selection q

	Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(n, n));
	thinQ = q * thinQ;
	end_q = omp_get_wtime();
	delta_q = end_q - start_q;
	cout << "Time for eigen version of Q: " << delta_q << endl;
	//cout << "Matrix q " << endl << thinQ << endl;
	//cout << "Matrix r " << endl;
	//cout << r << endl;

	check_with_eigen_result(qrAlgo->get_Q(), qrAlgo->get_R(), thinQ, r, n);
	// точность eigen в L2
	return 0;
}