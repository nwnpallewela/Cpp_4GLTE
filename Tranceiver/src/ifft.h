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

class IFFT {
private:
	double **ifftdata;
	fftw_complex in2[12], out[12];
							  //, c_out1[12], c_out2[12], c_in1[12], c_in2[12]; /* double [2] */
	fftw_plan q;
	//, cq1, cq2;

public:
	double** getifftw(int N, double **modulateddata, int len);
	double** getifftw_channel(int N, double **modulateddata, int len);
	IFFT(int len) {
		//cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);
		q = fftw_plan_dft_1d(12, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

		//cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 12; i++) {
			in2[i][0] = 0;
			in2[i][1] = 0;
			out[i][0] = 0;
			out[i][1] = 0;

	/*		c_out1[i][0] = 0;
			c_out1[i][1] = 0;
			c_out2[i][0] = 0;
			c_out2[i][1] = 0;

			c_in1[i][0] = 0;
			c_in1[i][1] = 0;
			c_in2[i][0] = 0;
			c_in2[i][1] = 0;*/
		}
		ifftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			ifftdata[i] = new double[len * 3 / 2];

			for (int j = 0; j < len; j++) {
				ifftdata[i][j] = 0;
			}
		}
	}
	IFFT() {
		int len = 24;
	//	cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);
		q = fftw_plan_dft_1d(12, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

	//	cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 12; i++) {
			in2[i][0] = 0;
			in2[i][1] = 0;
			out[i][0] = 0;
			out[i][1] = 0;

		/*	c_out1[i][0] = 0;
			c_out1[i][1] = 0;
			c_out2[i][0] = 0;
			c_out2[i][1] = 0;

			c_in1[i][0] = 0;
			c_in1[i][1] = 0;
			c_in2[i][0] = 0;
			c_in2[i][1] = 0;*/
		}
		ifftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			ifftdata[i] = new double[len * 3 / 2];

			for (int j = 0; j < len; j++) {
				ifftdata[i][j] = 0;
			}
		}
	}
};
double** IFFT::getifftw(int N, double **fftdata, int len) {
	/* double [2] */

	int count_ifftdata = 0;
	for (int j = 0; j < len; j = j + 8) {
		out[0][0] = 0;
		out[1][0] = 0;
		out[10][0] = 0;
		out[11][0] = 0;
		out[0][1] = 0;
		out[1][1] = 0;
		out[10][1] = 0;
		out[11][1] = 0;

		for (int i = 0; i < 8; i++) {
			out[i + 2][0] = fftdata[0][j + i];
			out[i + 2][1] = fftdata[1][j + i];
		}

		fftw_execute(q);

		/* normalize */
		for (int i = 0; i < N; i++) {

		//	cout<<in2[i][0]<<"*****"<<endl;
			in2[i][0] *= 1. / N;
			in2[i][1] *= 1. / N;

			if ((len * 3 / 2) > count_ifftdata) {
				ifftdata[0][count_ifftdata] = in2[i][0];
				ifftdata[1][count_ifftdata] = in2[i][1];
				count_ifftdata++;
			}

		}

		//		for (int i = 0; i < 12; i++){
		//				printf("recover: %3d %+9.5f %+9.5f I vs. %+9.5f %+9.5f I\n", i,
		//						out[i][0], out[i][1], in2[i][0], in2[i][1]);
		//			}

	}
	//fftw_destroy_plan(q);
	//fftw_cleanup();

	return ifftdata;
}


