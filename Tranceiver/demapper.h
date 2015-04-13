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

class DeMapper {

private:
	double** fft(double** tx);
	//double** ifft(double** tx);
	double **rxtemp1;
	double **rxtemp2;
	double **rx1;
	double **rx2;
	double *rxout;
	fftw_complex in[12], out[12]; /* double [2] */
	fftw_complex in2[8], out2[8];
	//fftw_plan q;
	fftw_plan pch;

public:
	double* get_demapped_rx(double** rx);
	DeMapper() {
		//	q = fftw_plan_dft_1d(8, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		pch = fftw_plan_dft_1d(12, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		rxtemp1 = new double*[2];
		rxtemp2 = new double*[2];
		rx1 = new double*[2];
		rx2 = new double*[2];
		rxout = new double[32];
		for (int i = 0; i < 2; i++) {
			rxtemp1[i] = new double[12];
			rxtemp2[i] = new double[12];
			rx1[i] = new double[12];
			rx2[i] = new double[12];

			for (int j = 0; j < 12; j++) {
				rxtemp1[i][j] = 0;
				rxtemp2[i][j] = 0;
				rx1[i][j] = 0;
				rx2[i][j] = 0;
			}
		}

		for (int i = 0; i < 12; i++) {
			in[i][0] = 0;
			in[i][1] = 0;
			out[i][0] = 0;
			out[i][1] = 0;

		}
		for (int i = 0; i < 8; i++) {
			in2[i][0] = 0;
			in2[i][1] = 0;
			out2[i][0] = 0;
			out2[i][1] = 0;

		}

	}

};

double* DeMapper::get_demapped_rx(double** rx) {

	for (int i = 0; i < 12; i++) {
		rx1[0][i] = rx[0][i];
		rx1[1][i] = rx[1][i];
		rx2[0][i] = rx[0][i + 12];
		rx2[1][i] = rx[1][i + 12];
	}

	/*double **rx1_t;
	 rx1_t=rx1;
	 rx1_t = fft(rx1);*/
	rxtemp1 = fft(rx1);
	rxtemp2 = fft(rx2);

	/*	double rxtemp1[][] = ifft(fft(rx1));// use if you are not using Equalizer in testing
	 double rxtemp2[][] = ifft(fft(rx2));
	 */

	for (int i = 0; i < 8; i++) {
		rxout[i] = rxtemp1[0][i + 2];
		rxout[i + 8] = rxtemp2[0][i + 2];
		rxout[i + 16] = rxtemp1[1][i + 2];
		rxout[i + 24] = rxtemp2[1][i + 2];

	}

	//////////////////////////////////////////////////////////////////use only you dont use Equalizer

	/*for (int i = 0; i < rxout.length; i++) {
	 double temp=(2*Math.round((Math.abs(rxout[i])-1)/2)+1);

	 if(Double.compare(temp, 7.0)>0){
	 temp =7.0;
	 }
	 if(Double.compare(rxout[i], 0.0)>0){

	 rxout[i]=temp;
	 }else{
	 rxout[i]= -temp;
	 }
	 }*/
	/////////////////////////////////////only run above code if you dont use equalizer
	////////////////////////////////////
	/*
	 * System.out.println("******************* rxout "); for (int i = 0; i <
	 * rxout.length; i++) { System.out.print(" "+i+" : "+rxout[i]+"\t"); }
	 */
	return rxout;
}

double** DeMapper::fft(double** tx) {

	for (int i = 0; i < 12; i++) {
		in[i][0] = tx[0][i];
		in[i][1] = tx[1][i];
	}

	/* forward Fourier transform, save the result in 'out' */
//	pch = fftw_plan_dft_1d(12, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(pch);

	for (int i = 0; i < 12; i++) {
		tx[0][i] = out[i][0];
		tx[1][i] = out[i][1];
	}

	return tx;
}
