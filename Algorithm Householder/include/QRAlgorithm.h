#pragma once

#include <omp.h>

namespace QR {

	void print_matrix(double* MatrixA, int Size);
	void multiplication(double* MatrixA, double* MatrixB, double* &MatrixC, int n);
class QRAlgorithm {
protected:
	int n;
	double* Q, *R;
	double* P; // matrix of v vectors;
	virtual double Norma(double* x, int size);

	virtual void copy_matrix(double* &destiny, double* &source);
	double* CreateMatrix_V(double* v);
public:
	double* A;
	QRAlgorithm(int _n);
	~QRAlgorithm()
	{
		delete[] A;
		delete[] Q;
		delete[] R;
		delete[] P;
	}
	void generator();
	virtual void QRDecomposition() = 0;
	virtual void QSelector() = 0;
	void check_result();
	double* get_Q() { return Q; }
	double* get_R() { return R; }
};

class PrimitiveQR : public QRAlgorithm{
public:
	PrimitiveQR(int _n) : QRAlgorithm(_n){};
	virtual void QRDecomposition();
	virtual void QSelector();
};

class RowHouseQR : public PrimitiveQR {
protected:
	void row_house(double*& _A, double* v, int size, int k);
public:
	RowHouseQR(int _n) : PrimitiveQR(_n) {};
	virtual void QRDecomposition();
};

class MultiplicationQR : public RowHouseQR {
public:
	MultiplicationQR(int _n) : RowHouseQR(_n) {};
	//void QRDecomposition();
	virtual void QSelector() override;
};

class FinalSequenceQR : public RowHouseQR {
public:
	FinalSequenceQR(int _n) : RowHouseQR(_n) {};
	void QRDecomposition();
    void QSelector() override;
};


class ParallelQR : public QRAlgorithm {
protected:
	virtual double Norma(double* x, int size) override;

	virtual void copy_matrix(double* &destiny, double* &source) override;
	void row_house(double*& _A, double* v, int size, int k);
public:
	ParallelQR(int _n) : QRAlgorithm(_n) {};

	void QRDecomposition() override;
	void QSelector() override;
};

class GivensRotation : public QRAlgorithm {
protected:
	void givens(double a, double b, double &c, double &s);
	void rowRotation(double*& A, double c, double s, int _i);
	void columnRotation(double*& A, double c, double s, int _j);
	double calculateP(double c, double s);
	void calculateCSfromP(double p, double &c, double &s);

	double *p;

public:
	GivensRotation(int _n) : QRAlgorithm(_n) {
		int size = (n*n - n) / 2;
		p = new double[size];
	};
	virtual void QRDecomposition();
	virtual void QSelector();
};

}