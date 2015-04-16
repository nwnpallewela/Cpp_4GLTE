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



using namespace std;
using namespace arma;



typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;


int findindex(string code);


class Modulator {

private:
	double **modulateddata;

public:
	double** modulate64QAM(string encodeddata);
	Modulator(int len) {
		modulateddata = new double*[2];
		for (int i = 0; i < 2; i++) {
			modulateddata[i] = new double[len];	//encodeddata.length() / 6
			for (int j = 0; j < len; j++) {
				modulateddata[i][j] = 0;
			}
		}
	}

};

double** Modulator::modulate64QAM(string encodeddata) {

//////////////////////////////////////////////////////////////////////////////////////////////////////////

	/*	double **modulateddata = new double*[2];
	 for (int i = 0; i < 2; i++) {
	 modulateddata[i] = new double[encodeddata.length() / 6];
	 for (int j = 0; j < encodeddata.length() / 6; j++) {
	 modulateddata[i][j] = 0;
	 }
	 }*/
	const int real[64] = { -7, -7, -7, -7, -7, -7, -7, -7, -5, -5, -5, -5, -5,
			-5, -5, -5, -1, -1, -1, -1, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3,
			-3, -3, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1,
			1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3 };
	;

	const int complex[64] = { -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5,
			1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7,
			-5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3,
			7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3 };

	int index;
	int count = 0;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j <= (encodeddata.length() - 6); j = j + 6) {

		index = findindex(encodeddata.substr(j, 6));
//cout<<real[index]<<" + "<<complex[index]<<endl;
		modulateddata[0][count] = real[index];
		modulateddata[1][count] = complex[index];
		count++;
//	cout<<encodeddata.substr(j, 6)<<" "<<real[index]<<" "<<complex[index]<<endl;
	}
	//delete[] real;
	//delete[] complex;
//cout<<endl<<"-------------------------"<<endl;
	return modulateddata;
}

int findindex(string code) {

	return ((code[0] - 48) + (code[1] - 48) * 2 + (code[2] - 48) * 4
			+ (code[3] - 48) * 8 + (code[4] - 48) * 16 + (code[5] - 48) * 32);
}
