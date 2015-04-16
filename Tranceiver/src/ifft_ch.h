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

class IFFT_CH {

private:
	double **ifftdata;
	fftw_complex c_out1[12],c_in1[12];;//, c_out2[12],  c_in2[12]; /* double [2] */
	fftw_plan cq1;//, cq2;

public:
	//double** getifftw(int N, double **modulateddata, int len);
	double** getifftw_channel(int N, double **modulateddata, int len);
	IFFT_CH(int len) {
		cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);

		//cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		for (int i = 0; i < 12; i++) {

			c_out1[i][0] = 0;
			c_out1[i][1] = 0;
		/*	c_out2[i][0] = 0;
			c_out2[i][1] = 0;*/

			c_in1[i][0] = 0;
			c_in1[i][1] = 0;
			/*c_in2[i][0] = 0;
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
	IFFT_CH() {
		int len = 24;

		for (int i = 0; i < 12; i++) {

			c_out1[i][0] = 0;
			c_out1[i][1] = 0;
			/*c_out2[i][0] = 0;
			c_out2[i][1] = 0;*/

			c_in1[i][0] = 0;
			c_in1[i][1] = 0;
		/*	c_in2[i][0] = 0;
			c_in2[i][1] = 0;*/
		}
		ifftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			ifftdata[i] = new double[len * 3 / 2];

			for (int j = 0; j < len; j++) {
				ifftdata[i][j] = 0;
			}
		}
		cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);
	//	cq2 = fftw_plan_dft_1d(12, c_out2, c_in2, FFTW_BACKWARD, FFTW_ESTIMATE);
	}
};

double** IFFT_CH::getifftw_channel(int N, double **fftdata, int len) {
	/* double [2] */
	cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);

	int count_ifftdata = 0;
	for (int j = 0; j < 2; ++j) {

		for (int i = 0; i < 12; i++) {
			c_out1[i][0] = fftdata[0][i + 12 * j];
			c_out1[i][1] = fftdata[1][i + 12 * j];
			/*c_out2[i][0] = fftdata[0][i + 12];
			 c_out2[i][1] = fftdata[1][i + 12];*/

			//std::cout << "" << c_out1[i][0] << std::endl << c_out1[i][1] << std::endl;
		}
//	std::cout << "******************" << std::endl;
		//fftw_destroy_plan(cq1);
		//fftw_cleanup();

//	fftw_destroy_plan(cq1);
//	fftw_cleanup();

		fftw_execute(cq1);
		for (int i = 0; i < N; i++) {
			//	std::cout << c_in1[i][0] << std::endl << c_in1[i][1] << std::endl;
			//	std::cout << c_in2[i][0] << std::endl << c_in2[i][1] << std::endl;
			c_in1[i][0] *= 1. / N;
			c_in1[i][1] *= 1. / N;

			if ((len * 3 / 2) > count_ifftdata) {

				ifftdata[0][count_ifftdata] = c_in1[i][0];
				ifftdata[1][count_ifftdata] = c_in1[i][1];
				count_ifftdata++;
			}

		}
	}

	/*	count_ifftdata = 0;

	 cq1 = fftw_plan_dft_1d(12, c_out1, c_in1, FFTW_BACKWARD, FFTW_ESTIMATE);

	 fftw_execute(cq1);

	 normalize
	 for (int i = 0; i < N; i++) {
	 //	std::cout << c_in1[i][0] << std::endl << c_in1[i][1] << std::endl;
	 //	std::cout << c_in2[i][0] << std::endl << c_in2[i][1] << std::endl;

	 c_in1[i][0] *= 1. / N;
	 c_in1[i][1] *= 1. / N;

	 if ((len * 3 / 2) > count_ifftdata) {
	 ifftdata[0][count_ifftdata] = c_in1[i][0];
	 ifftdata[1][count_ifftdata] = c_in1[i][1];

	 count_ifftdata++;
	 }

	 }*/

	//fftw_destroy_plan(cq1);
	//	fftw_cleanup();
	//		for (int i = 0; i < 12; i++){
	//				printf("recover: %3d %+9.5f %+9.5f I vs. %+9.5f %+9.5f I\n", i,
	//						out[i][0], out[i][1], in2[i][0], in2[i][1]);
	//			}
	//fftw_destroy_plan(q);
	//fftw_cleanup();
	/*std::cout << "                     " << ifftdata[0][0] << " " << c_in1[0][0]
	 << std::endl;*/
	return ifftdata;
}
