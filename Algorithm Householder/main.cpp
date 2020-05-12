#include <iostream>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>
#include "include/QRAlgorithm.h"

using namespace std;
using namespace QR;
//using namespace Eigen;

void check_with_eigen_result(double* MatrixQ, double* MatrixR, Eigen::MatrixXd MatrixQE, Eigen::MatrixXd MatrixRE, int size) {
	cout << "Check result matrices with eigen results... " << endl;
	bool flag_ok = true;
	cout << "Check matrix R: " << endl;
	double eps = pow(10, -7); //10 ^ (-7);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i < j)
				if (abs(MatrixR[i* size + j] - MatrixRE(i, j)) > eps) {
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
			if (abs(MatrixQ[i* size + j] - MatrixQE(i, j)) > eps) {
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

//int main() {
int main(int argc, char **argv) {
	//cout << "Please choose mode: 0 - Test mode with definitly 3x3 matrix, 1 - primitive version, 2 - rowHouse version, 3 - multipilication version, 4 - final sequence version" << endl;

	int n, mode;
	mode = atoi(argv[1]);
	n = atoi(argv[2]);

	//cin >> mode;
	//cout << "Please enter size of matrix: ";
	//cin >> n;
	//cout << endl;

	//mode = 7;
	//n = 1000;
	cout.scientific;
	QRAlgorithm *qrAlgo;

	switch (mode) {
	case -1:
	{
		// example from wiki
		qrAlgo = new GivensRotation(n);
		qrAlgo->A[0] = 6;	qrAlgo->A[1] = 5;	qrAlgo->A[2] = 0;
		qrAlgo->A[3] = 5; 	qrAlgo->A[4] = 1;	qrAlgo->A[5] = 4;
		qrAlgo->A[6] = 0; 	qrAlgo->A[7] = 4; 	qrAlgo->A[8] = 3;
		break;
	}
	case 0:
	{
		qrAlgo = new GivensRotation(n);
		qrAlgo->A[0] = 1;	qrAlgo->A[1] = -2;	qrAlgo->A[2] = 1;
		qrAlgo->A[3] = 2; 	qrAlgo->A[4] = -1;	qrAlgo->A[5] = 3;
		qrAlgo->A[6] = 2; 	qrAlgo->A[7] = 3; 	qrAlgo->A[8] = -4;
		break;
	}
	case 1: 
	{
		qrAlgo = new PrimitiveQR(n);
		break;
	}
	case 2:
	{
		qrAlgo = new RowHouseQR(n);
		break;
	}
	case 3: 
	{
		qrAlgo = new MultiplicationQR(n);
		break;
	}
	case 4:
	{
		qrAlgo = new FinalSequenceQR(n);
		break;
	}
	case 5:
	{
		int count_threads = 2;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelQR(n);
		break;
	}
	case 6:
	{
		int count_threads = 4;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelQR(n);
		break;
	}
	case 7:
	{
		int count_threads = 8;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelQR(n);
		break;
	}
	case 8:
	{
		qrAlgo = new GivensRotation(n);
		break;
	}
	case 9:
	{
		int count_threads = 2;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelGivensRotation(n);
		break;
	}
	case 10:
	{
		int count_threads = 4;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelGivensRotation(n);
		break;
	}
	case 11:
	{
		int count_threads = 8;
		omp_set_num_threads(count_threads);
		qrAlgo = new ParallelGivensRotation(n);
		break;
	}
	default:
		break;

	}

	qrAlgo->generator();

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

	//Eigen::MatrixXd m(n, n);
	//for (int i = 0; i < n; i++) {
	//	for (int j = 0; j < n; j++)
	//		m(i, j) = qrAlgo->A[i*n + j];
	//}
	//auto qr = m.householderQr(); // casted type Matrix -> Householder matrix

	//start_decomp = omp_get_wtime();
	//
	//auto r = qr.matrixQR(); // qr decomposition, r - upper triangle, vectors v - lower triangle

	//end_decomp = omp_get_wtime();
	//delta_decomp = end_decomp - start_decomp;
	//cout << "Time for eigen version of decomposition: " << delta_decomp << endl;

	//start_q = omp_get_wtime();
	//auto q = qr.householderQ(); // selection q

	//Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(n, n));
	//thinQ = q * thinQ;
	//end_q = omp_get_wtime();
	//delta_q = end_q - start_q;
	//cout << "Time for eigen version of Q: " << delta_q << endl;
	//cout << "Matrix q " << endl << thinQ << endl;
	//cout << "Matrix r " << endl;
	//cout << r << endl;

	//check_with_eigen_result(qrAlgo->get_Q(), qrAlgo->get_R(), thinQ, r, n);
	// точность eigen в L2
	return 0;
}