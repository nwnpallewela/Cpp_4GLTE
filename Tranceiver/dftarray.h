#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>
#include <complex>
#include <valarray>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <cstdlib>
#include <armadillo>

const double PI = 3.141592653589793238460;

class DFTARRAY {

private:
	arma::cx_mat DF;

public:

	arma::cx_mat getDFT();
	void printDFT();
	DFTARRAY() {
		DF = zeros < cx_mat > (16, 16);

		arma::cx_double val = cx_double(0.0, 0.0);
		for (int i = 0; i < 8; i++) {	//df.size()/2
			for (int j = 0; j < 8; j++) {
				val = exp(cx_double(0.0,(-1 * (2 * PI * (i) * (j) / 8))));

				if (abs(val.real()) < pow(10, -10)) {
					val = cx_double(0.0, val.imag());
				}
				if (abs(val.imag()) < pow(10, -10)) {
					val = cx_double(val.real(), 0.0);
				}
				DF.at(i, j) = val;

				DF.at(i + 8, j + 8) = val;
				DF.at(i, j + 8) = cx_double(0.0, 0.0);
				DF.at(i + 8, j) = cx_double(0.0, 0.0);
			}
		}

	}

};

arma::cx_mat DFTARRAY::getDFT() {
	return DF;
}

void DFTARRAY::printDFT() {
	cout << "DFT Matrix :" << endl;
	for (int i = 0; i < DF.size(); i++) { //16
		for (int j = 0; j < DF.size(); j++) {
			cout << DF(i, j) << "    ";
		}
		cout << endl;
	}

}
