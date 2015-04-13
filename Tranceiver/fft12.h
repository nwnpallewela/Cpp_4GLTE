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

class FFT_12 {

private:
	double **fftdata;
//	fftw_complex inch[12], outch[12]; /* double [2] */
	fftw_complex *inch, *outch;
	fftw_plan p2;

public:
	double** getfftw(int N, double **txdata, int len);

	FFT_12() {
		inch = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 12);
		outch = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 12);

		fftw_complex temp[12];
		p2 = fftw_plan_dft_1d(12, inch, outch, FFTW_FORWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 12; i++) {
			inch[i][0] = 0;
			inch[i][1] = 0;
			outch[i][0] = 0;
			outch[i][1] = 0;
		}
		fftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			fftdata[i] = new double[12];

			for (int j = 0; j < 12; j++) {
				fftdata[i][j] = 0;
			}
		}
	}

};

double** FFT_12::getfftw(int N, double **txdata, int len) {


	/* forward Fourier transform, save the result in 'out' */
	//p2 = fftw_plan_dft_1d(12, inch, outch, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < 12; i++) {
			inch[i][0] = 0;
			inch[i][1] = 0;
			outch[i][0] = 0;
			outch[i][1] = 0;
		}
	for (int i = 0; i < N; i++) {
		inch[i][0] = txdata[0][i];
		inch[i][1] = txdata[1][i];

	}
	fftw_execute(p2);

	for (int i = 0; i < N; i++) {
		fftdata[0][i] = outch[i][0];
		fftdata[1][i] = outch[i][1];
		//cout<<out[i][0]<<"...";
	}

	// fftw_destroy_plan(p2);
//fftw_cleanup();

	return fftdata;

}
