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

class FFT {

private:
	double **fftdata;
	fftw_complex in[8], out[8]; /* double [2] */

	fftw_plan p;

public:
	double** getfftw(int N, double **modulateddata, int len);
	FFT(int len) {
		fftw_complex temp[8];
		p = fftw_plan_dft_1d(8, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 8; i++) {
			in[i][0] = 0;
			in[i][1] = 0;
			out[i][0] = 0;
			out[i][1] = 0;
		}
		fftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			fftdata[i] = new double[len];

			for (int j = 0; j < len; j++) {
				fftdata[i][j] = 0;
			}
		}
	}

};

double** FFT::getfftw(int N, double **modulateddata, int len) {

	for (int j = 0; j < len; j = j + 8) {

		for (int i = 0; i < N; i++) {
			in[i][0] = modulateddata[0][j + i];
			in[i][1] = modulateddata[1][j + i];
		}

		/* forward Fourier transform, save the result in 'out' */

		fftw_execute(p);
		for (int i = 0; i < N; i++) {
			fftdata[0][i + j] = out[i][0];
			fftdata[1][i + j] = out[i][1];
		}

	}
	//fftw_destroy_plan(p);

	return fftdata;

}
