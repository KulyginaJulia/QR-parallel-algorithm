#include <iostream>
#include <list>
#include <vector>
#include <cmath>

using namespace std;

double Normirovka(vector<double> s, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += s[i] * s[i];
	return sum;
}

vector<vector<double>> CreateMatrix_V(vector<double> v) {
	int size = v.size();
	vector<vector<double>> V(size);
	double norm = Normirovka(v, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double val = v[i] * v[j] / norm;
			V[i].push_back(val);
		}
	}

	return V;
}

vector<vector<double>> multiplication(vector<vector<double>> MatrixA, vector<vector<double>> MatrixB) {
	int N = MatrixA.size();
	vector<vector<double>>MatrixC(N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrixC[j].push_back(0);
			for (int k = 0; k < N; k++) {
				double val = MatrixC[j][i] + MatrixA[i][k] * MatrixB[j][k];
				MatrixC[j][i] = val;
			}
		}
	}
	return MatrixC;
}

int main() {
	vector<vector<double>> C{ { 1, 2, 2 },{ -2, -1, 3 },{ 1, 3, -4 } };
	vector<vector<double>> E{ { 1, 0, 0 },{ 0, 1, 0 },{ 0, 0, 1 } };

	int n = 3;

	list< vector<vector<double>>> list_P;
	for (int k = 0; k < n - 1; k++) {
		vector<double> v;
		for (int i = 0; i < n; i++) {
			int signum = C[k][k] > 0 ? 1 : -1;
			if (i < k) {
				v.push_back(0);
			}
			else {
				v.push_back(C[k][i] + signum * sqrt(Normirovka(C[k], n)) * E[k][i]);
			}

		}
		vector<vector<double>> V = CreateMatrix_V(v);
		vector<vector<double>> P1(n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double val = E[i][j] - 2 * V[i][j];
				P1[i].push_back(val);
			}
		}
		list_P.push_back(P1);
		C = multiplication(P1, C);
		P1.clear();
		v.clear();
	}
	vector<vector<double>> Q(n);
	Q = multiplication(list_P.back(), list_P.front());

	return 0;
}