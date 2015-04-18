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

class IFFT_Z {
private:
	double **ifftdata;
	fftw_complex in2[8], out[8];
							  //, c_out1[12], c_out2[12], c_in1[12], c_in2[12]; /* double [2] */
	fftw_plan q;
	//, cq1, cq2;

public:
	double** getifftw(int N, arma::mat Y, int len);
	double** getifftw_channel(int N, double **modulateddata, int len);
	IFFT_Z(int len) {
		//cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);
		q = fftw_plan_dft_1d(8, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

		//cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 8; i++) {
			in2[i][0] = 0;
			in2[i][1] = 0;
			out[i][0] = 0;
			out[i][1] = 0;


		}
		ifftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			ifftdata[i] = new double[len];

			for (int j = 0; j < len; j++) {
				ifftdata[i][j] = 0;
			}
		}
	}
	IFFT_Z() {
		int len = 16;
	//	cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);
		q = fftw_plan_dft_1d(8, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

	//	cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 8; i++) {
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
			ifftdata[i] = new double[len ];

			for (int j = 0; j < len; j++) {
				ifftdata[i][j] = 0;
			}
		}
	}
};
double** IFFT_Z::getifftw(int N, arma::mat Y, int len) {
	/* double [2] */

	int count_ifftdata = 0;
	for (int j = 0; j < len; j = j + 8) {
	//	cout<<"This is iift in"<<endl;
		for (int i = 0; i < 8; i++) {
			out[i][0] = Y.at(j+i,0);//[0][j + i];
			out[i][1] =Y.at(j+i+16,0);
		//	cout<<out[i][0]<<" + "<<out[i][1]<<endl;
		}

		fftw_execute(q);
	//	cout<<"This is iift out"<<endl;
		/* normalize */
		for (int i = 0; i < N; i++) {

		//	cout<<in2[i][0]<<"*****"<<endl;
			in2[i][0] *= 1. / N;
			in2[i][1] *= 1. / N;

		//	cout<<in2[i][0]<<" + "<<in2[i][1]<<endl;
				ifftdata[0][count_ifftdata] = in2[i][0];
				ifftdata[1][count_ifftdata] = in2[i][1];
				count_ifftdata++;


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


