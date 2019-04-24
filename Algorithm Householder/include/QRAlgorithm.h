#pragma once



class QRAlgorithm {
protected:
	int n;
	double** Q, **R;
	double** P; // matrix of v vectors;
	double Norma(double* x);
	double** multiplication(double** MatrixA, double** MatrixB);
	void copy_matrix(double** &destiny, double** &source);
	double** CreateMatrix_V(double* v);
public:
	double** A;
	QRAlgorithm(int _n);
	void generator();
	virtual void QRDecomposition() = 0;
	virtual void QSelector() = 0;
	void check_result();
	double** get_Q() { return Q; }
	double** get_R() { return R; }
};

class PrimitiveQR : public QRAlgorithm{
public:
	PrimitiveQR(int _n) : QRAlgorithm(_n){};
	void QRDecomposition();
	void QSelector();
};

class RowHouseQR : public PrimitiveQR {
protected:
	double** row_house(double** _A, double* v);
public:
	RowHouseQR(int _n) : PrimitiveQR(_n) {};
	void QRDecomposition();
};
//
//class MultiplicationQR : QRAlgorithm {
//public:
//	void QRDecomposition();
//	void QSelector();
//};