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

#include "ifft_ch.h"

class Channel {
private:
	//Complex tr;
	string channel_type;
	string corr_type;
	double fcarry;
	double sigma;
	double tx_corr_coeff, rx_corr_coeff;
	int dopp_freq;
	int no_taps;
	arma::mat sqrt_corr_matrix;
	arma::mat tx_corr_matrix;
	arma::mat rx_corr_matrix;
	arma::mat corr_matrix;
	arma::cx_mat H;
	arma::cx_mat H_r;
	arma::cx_double tr_1_coeff;
	arma::cx_double tr_2_coeff;
	arma::cx_double tr_1_calc;
	arma::cx_double tr_2_calc;
	arma::cx_mat A;
	arma::cx_mat B;
	arma::cx_mat noise;
//	FFT_12 fft;
	//IFFT_CH ifft;
	arma::cx_mat R;
	arma::cx_mat Tx;
	double **Rx;
	arma::cx_mat H1_mat;
	arma::cx_mat DF;
	arma::cx_mat H_out;
	DFTARRAY dft;

public:
	double** LTEMIMOCHANNEL(double **Tx1, int tx1length, double **Tx2);
	cx_mat LTEMIMOCHANNEL_mat();
	arma::cx_mat getHout();
	Channel(string corr_type, string channel_type, double fcarry,
			double sigma) {
		this->channel_type = channel_type;
		this->corr_type = corr_type;
		this->fcarry = fcarry;
		this->sigma = sigma;
		this->dopp_freq = 5;
		this->no_taps = 7;
		this->tx_corr_coeff = 0;
		this->rx_corr_coeff = 0;
		sqrt_corr_matrix = mat(4, 4, fill::zeros);
		tx_corr_matrix = mat(2, 2, fill::zeros); //correlation matrix initiallization
		rx_corr_matrix = mat(2, 2, fill::zeros);
		corr_matrix = kron(tx_corr_matrix, rx_corr_matrix);
		H = zeros < cx_mat > (24, 24);
		H_r= zeros < cx_mat > (24, 24);
		H.at(0,0)=cx_double(4.292090270076418,-17.645757923301453);
		H.at(0,12)=cx_double(5.7065576483675144,4.956916362581227);
		H.at(12,0)=cx_double(8.462940503301988,-8.556251748934612);
		H.at( 12,12)=cx_double(17.92479447394951,3.726255843428431);

		H.at(1,1)=cx_double(-11.984837606866979,8.730609470686591);
		H.at(1,13)=cx_double(-0.029494809647753306,2.706869698671653);
		H.at(13,1)=cx_double(11.522949390771974,-2.7406697799707116);
		H.at( 13,13)=cx_double(-10.86972574330876,-12.300258015190952);

		H.at(2,2)=cx_double(-6.334994589122318,4.8896610723011635);
		H.at(2,14)=cx_double(-0.8419481818122718,3.4639762532035943);
		H.at(14,2)=cx_double(5.312419919882708,5.627346467335482);
		H.at( 14,14)=cx_double(-7.720635676164648,1.3003711382188823);

		H.at(3,3)=cx_double(-5.919776236225726,-4.138949525054658);
		H.at(3,15)=cx_double(-4.2320103591584,-1.6673722136773326);
		H.at(15,3)=cx_double(4.019677449652807,-9.877341163794858);
		H.at( 15,15)=cx_double(-6.564524935820653,7.638217069039449);

		H.at(4,4)=cx_double(12.677438392041097,-11.91270168461972);
		H.at(4,16)=cx_double(13.789453371572277,-1.8414492065988721);
		H.at(16,4)=cx_double(-16.831847689682583,-11.4907075075654);
		H.at( 16,16)=cx_double(-12.804285858885985,1.9946656974420858);

		H.at(5,5)=cx_double(10.739119494190227,-10.59260050546974);
		H.at(5,17)=cx_double(2.135394146412287,-7.634413387467054);
		H.at(17,5)=cx_double(13.151930364984523,10.36632617096652);
		H.at( 17,17)=cx_double(-0.5004646047298407,5.471507799446753);

		H.at(6,6)=cx_double(11.171387740013202,-5.676360518103919);
		H.at(6,18)=cx_double(-3.186610937958143,-5.594051911145665);
		H.at(18,6)=cx_double(-3.8669736555478575,-13.086482906575613);
		H.at( 18,18)=cx_double(2.762515381340146,1.9215139260204825);

		H.at(7,7)=cx_double(8.596220190684601,9.38896241612708);
		H.at(7,19)=cx_double(-3.5805947881882214,-18.568765419081487);
		H.at(19,7)=cx_double(5.903952652466867,2.300952042701443);
		H.at( 19,19)=cx_double(1.8751860738080768,-2.1231543003644253);

		H.at(8,8)=cx_double(-17.284536383505248,-12.852770831628764);
		H.at(8,20)=cx_double(-7.7609521231337375,-11.208697817753587);
		H.at(20,8)=cx_double(-12.905433972189195,-0.11404978698001234);
		H.at( 20,20)=cx_double(-4.7382763176434946,13.205777526684454);

		H.at(9,9)=cx_double(-12.722325418495084,-1.8919863492961282);
		H.at(9,21)=cx_double(3.9511518755150146,0.9638594385724886);
		H.at(21,9)=cx_double(-1.9788354778521713,14.530435805648807);
		H.at( 21,21)=cx_double(3.836863319170617,10.991156548155512);

		H.at(10,10)=cx_double(-5.18859734851744,6.824672406885832);
		H.at(10,22)=cx_double(-16.15246188666331,4.481227906838198);
		H.at(22,10)=cx_double(-7.112660071246864,18.188903294637505);
		H.at( 22,22)=cx_double(-11.546091571754225,-22.187265885837707);

		H.at(11,11)=cx_double(-7.3947624958475515,-0.11400181978302149);
		H.at(11,23)=cx_double(-4.084706318034912,-5.775298736083231);
		H.at(23,11)=cx_double(4.651255406879632,2.9880385433609415);
		H.at( 23,23)=cx_double(3.326598589706287,-5.963794114434786);
		tr_1_coeff = cx_double(1.0, 0.0);
		tr_2_coeff = cx_double(1.0, 0.0);
		tr_1_calc = cx_double(0.0, 0.0);
		tr_2_calc = cx_double(0.0, 0.0);
		A = randn < cx_mat > (4, 1);
		B = sqrt_corr_matrix * A;
		noise = randn < cx_mat > (24, 1);
//		fft = FFT_12();
//		ifft = IFFT_CH();
		Tx = zeros < cx_mat > (24, 1);
		H1_mat = zeros < cx_mat > (16, 16);
		dft = DFTARRAY();

		DF = dft.getDFT();

		Rx = new double*[2];
		for (int i = 0; i < 2; i++) {
			Rx[i] = new double[24];

			for (int j = 0; j < 24; j++) {
				Rx[i][j] = 0;
			}
		}
	}

};

/*double** Channel::ifft(double **Tx1) {


 }*/
arma::cx_mat Channel::getHout() {
	return H_out;
}
cx_mat Channel::LTEMIMOCHANNEL_mat(){
	double Medium_type[4][4] = { { 0.84, 0.52, 0.13, 0.08 }, { 0.52, 0.84, 0.08,
				0.13 }, { 0.13, 0.08, 0.84, 0.52 }, { 0.08, 0.13, 0.52, 0.84 } };

		double High_type[4][4] = { { 0.72, 0.45, 0.45, 0.28 }, { 0.45, 0.72, 0.28,
				0.45 }, { 0.45, 0.28, 0.72, 0.45 }, { 0.28, 0.45, 0.45, 0.72 } };

		sqrt_corr_matrix = mat(4, 4, fill::zeros);
		if (corr_type == ("Low")) {
			tx_corr_coeff = 0;
			rx_corr_coeff = 0;
			for (int i = 0; i < 4; ++i) {
				sqrt_corr_matrix(i, i) = 1;
			}
		} else if (corr_type == ("Medium")) {
			tx_corr_coeff = 0.3;
			rx_corr_coeff = 0.9;
			for (int j = 0; j < 4; ++j) {
				for (int i = 0; i < 4; ++i) {
					sqrt_corr_matrix(j, i) = Medium_type[j][i];
				}
			}
		} else {
			tx_corr_coeff = 0.9;
			rx_corr_coeff = 0.9;
			for (int j = 0; j < 4; ++j) {
				for (int i = 0; i < 4; ++i) {
					sqrt_corr_matrix(j, i) = High_type[j][i];
				}
			}
		}
		tx_corr_matrix = mat(2, 2, fill::zeros); //correlation matrix initiallization
		rx_corr_matrix = mat(2, 2, fill::zeros);
		tx_corr_matrix(0, 0) = 1;
		tx_corr_matrix(0, 1) = tx_corr_coeff;
		tx_corr_matrix(1, 0) = tx_corr_coeff;
		tx_corr_matrix(1, 1) = 1;
		rx_corr_matrix(0, 0) = 1;
		rx_corr_matrix(0, 1) = rx_corr_coeff;
		rx_corr_matrix(1, 0) = rx_corr_coeff;
		rx_corr_matrix(1, 1) = 1;
		corr_matrix = kron(tx_corr_matrix, rx_corr_matrix);

		int path_delays[9] = { 0, 30, 70, 90, 110, 190, 410, 0, 0 };
		double path_gains[9] = { 0, -1, -2, -3, -8, -17.2, -20.8, 0, 0 };
		if (channel_type == ("EPA 5Hz")) {

			dopp_freq = 5;
			no_taps = 7;
		} else if (channel_type == ("EVA 5Hz")) {
			path_delays[2] = 150;
			path_delays[3] = 310;
			path_delays[4] = 370;
			path_delays[5] = 710;
			path_delays[6] = 1090;
			path_delays[7] = 1730;
			path_delays[8] = 2510;

			path_gains[1] = -1.5;
			path_gains[2] = -1.4;
			path_gains[3] = -3.6;
			path_gains[4] = -0.6;
			path_gains[5] = -9.1;
			path_gains[6] = -7;
			path_gains[7] = -12;
			path_gains[8] = -16.9;

			dopp_freq = 5;
			no_taps = 9;
		} else if (channel_type == ("EVA 70Hz")) {

			path_delays[2] = 150;
			path_delays[3] = 310;
			path_delays[4] = 370;
			path_delays[5] = 710;
			path_delays[6] = 1090;
			path_delays[7] = 1730;
			path_delays[8] = 2510;

			path_gains[1] = -1.5;
			path_gains[2] = -1.4;
			path_gains[3] = -3.6;
			path_gains[4] = -0.6;
			path_gains[5] = -9.1;
			path_gains[6] = -7;
			path_gains[7] = -12;
			path_gains[8] = -16.9;

			dopp_freq = 70;
			no_taps = 9;
		} else if (channel_type == ("ETU 70Hz")) {
			path_delays[1] = 50;
			path_delays[2] = 120;
			path_delays[3] = 200;
			path_delays[4] = 230;
			path_delays[5] = 500;
			path_delays[6] = 1600;
			path_delays[7] = 2300;
			path_delays[8] = 5000;

			path_gains[0] = -1;
			path_gains[1] = -1;
			path_gains[2] = -1;
			path_gains[3] = 0;
			path_gains[4] = 0;
			path_gains[5] = 0;
			path_gains[6] = -3;
			path_gains[7] = -5;
			path_gains[8] = -7;

			dopp_freq = 70;
			no_taps = 9;
		} else if (channel_type == ("ETU 300Hz")) {
			path_delays[1] = 50;
			path_delays[2] = 120;
			path_delays[3] = 200;
			path_delays[4] = 230;
			path_delays[5] = 500;
			path_delays[6] = 1600;
			path_delays[7] = 2300;
			path_delays[8] = 5000;

			path_gains[0] = -1;
			path_gains[1] = -1;
			path_gains[2] = -1;
			path_gains[3] = 0;
			path_gains[4] = 0;
			path_gains[5] = 0;
			path_gains[6] = -3;
			path_gains[7] = -5;
			path_gains[8] = -7;
			dopp_freq = 300;
			no_taps = 9;
		} else {
			cout << "ERROR : Unknown channel type " << endl;

			dopp_freq = 300;
			no_taps = 9;
		}

	///////////////////////////////////////////////////////////////////
		/*
		 cout << tx_corr_matrix.t() << endl;
		 cout << rx_corr_matrix.t() << endl;

		 cout << corr_matrix.t() << endl;
		 cout << sqrt_corr_matrix.t() << endl;
		 cout << (sqrt_corr_matrix * sqrt_corr_matrix).t() << endl;

		 cout << (sqrt_corr_matrix % sqrt_corr_matrix).t() << endl;
		 cout << (sqrt_corr_matrix * sqrt_corr_matrix).t() << endl;*/

	///////////////////////////////////////////////////////
		int l = 12;

		double f[l * 2];
		for (int k = 0; k < l; ++k) {

			f[k] = fcarry - 59 * 15 * 0.000001 + 15 * 0.000001 * k;
			f[l + k] = fcarry - 59 * 15 * 0.000001 + 15 * 0.000001 * k;

		}

		H_r = zeros < cx_mat > (2 * l, 2 * l); //zeros(2 * l) + 1i * zeros(2 * l);
		//cout << " Initial H matrix value :" << endl << H.t() << endl;
		//for k = 1:l
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*double tr_1_coeff = 1;						// can be complex...............
		 double tr_2_coeff = 1;*/
		//arma::cx_double test=
		tr_1_coeff = cx_double(1.0, 0.0);
		tr_2_coeff = cx_double(1.0, 0.0);
		tr_1_calc = cx_double(0.0, 0.0);
		tr_2_calc = cx_double(0.0, 0.0);
		A = randn < cx_mat > (4, 1);
		B = sqrt_corr_matrix * A;
		/*
		 cout << " Last A matrix value :" << endl << A.t() << endl;
		 cout << " Last B matrix value :" << endl << B.t() << endl;
		 */

		///////////////////////////////////////////////////////////////////////////////////////////////////////
		for (int k = 0; k < l; ++k) {
			A = randn < cx_mat > (4, 1);	//randn(4, 1) + 1i * randn(4, 1);
			B = sqrt_corr_matrix * A;

			tr_1_coeff = cx_double(1.0, 0.0);
			tr_2_coeff = cx_double(1.0, 0.0);
	//cout<<"EXP : "<<tr_1_coeff(0, 0)<<" : "<<exp(tr_1_coeff(0, 0))<<endl;
			// for m =1:no_taps
			for (int m = 0; m < no_taps; ++m) {
				tr_1_calc = cx_double(0.0, 2.0 * PI * f[k] * path_delays[m]);
				tr_2_calc = cx_double(0.0, 2.0 * PI * f[k + 8] * path_delays[m]);

				//	cout << "************** " << m << " " << k <<" EXP : "<<(tr_1_calc(0, 0))<<" : "<<exp(tr_1_calc(0, 0))<< endl;
				tr_1_coeff = tr_1_coeff
						+ sqrt(pow(10.0, (path_gains[m]))) * std::exp(tr_1_calc);
				tr_2_coeff = tr_2_coeff
						+ sqrt(pow(10.0, (path_gains[m]))) * std::exp(tr_2_calc);
			}
			//Complex.valueOf(10.0, 0).pow(path_gains[m]).sqrt()
			//.plus(tr_1_calc.exp()).plus(tr_1_coeff);

			H_r.at(k, k) = tr_1_coeff * B.at(0, 0);        //2 by 2 MIMO --> 4 Paths
			H_r.at(k, k + l) = tr_2_coeff * B.at(1, 0);
			H_r.at(k + l, k) = tr_1_coeff * B.at(2, 0);
			H_r.at(k + l, k + l) = tr_2_coeff * B.at(3, 0);

		}
		//H.print();
		/////////////////

		//////////////////
		/*cout << " tr_1_coeff(0, 0) * B(0, 0) : " << tr_1_coeff(0, 0) << " * "
		 << B(0, 0) << " = " << tr_1_coeff(0, 0) * B(0, 0) << endl;
		 cout << " Last A matrix value :" << endl << A.t() << endl;
		 cout << " Last B matrix value :" << endl << B.t() << endl;

		 cout << H.t() << endl;*/
		///////////////////////////////////////////////////////////////
	/*	noise = randn < cx_mat > (2 * l, 1);

		double **Tx_test;
		Tx_test = Tx1;
		//FFT_12 fft = FFT_12();
		Tx_test = fft.getfftw(12, Tx_test, 12);
		Tx_test = fft.getfftw(12, Tx1, 12);
		Tx1 = fft.getfftw(12, Tx1, 12);
		Tx2 = fft.getfftw(12, Tx2, 12);


		for (int var = 0; var < 12; ++var) {
			Tx.at(var, 0) = cx_double(Tx1[0][var], Tx1[1][var]);
			Tx.at(var + 12, 0) = cx_double(Tx2[0][var], Tx2[1][var]);
			//cout<<Tx1[0][var]<<"+"<< Tx1[1][var]<<"i * ";
		}
		Tx.print("This is channel tx: ");
		//cout<<"This is TX vector "<< Tx.t()<<endl;

	///////////////////////////////////////////////////////////////////////////////////////////
		R = H * Tx;	// + sigma * noise;

		//cout<<"This is R vector "<< R.t()<<endl;

		for (int var = 0; var < 24; ++var) {
			Rx[0][var] = R.at(var, 0).real();
			Rx[1][var] = R.at(var, 0).imag();

			//cout<<Rx[0][var]<<"+"<<Rx[1][var]<<"i \t";
		}
	//	cout << "this is RX b: " << Rx[0][0] << endl;
		Rx = ifft.getifftw_channel(12, Rx, 24);*/

	//	cout << "this is RX a: " << Rx[0][0] << endl;
		//Rx = ifft(Rx);
	//	arma::cx_mat Rx1 = ifft(R(1:8));
	//	arma::cx_mat Rx2 = ifft(R(9:16));
		///////////////////////////////////////////////////////////////////////////////////////

		//arma::cx_mat H_1 = zeros<cx_mat>(16, 16);

		//for y=1:8,
		for (int y = 0; y < 8; ++y) {
			for (int z = 0; z < 8; ++z) {
				//		H_1(y, z) = H(y + 2, z + 2);
				H1_mat.at(y, z) = H(y + 2, z + 2);
			}
		}

		for (int y = 0; y < 8; ++y) {
			for (int z = 8; z < 16; ++z) {
				H1_mat.at(y, z) = H(y + 2, z + 6);
			}
		}

		for (int y = 8; y < 16; ++y) {
			for (int z = 0; z < 8; ++z) {
				H1_mat.at(y, z) = H(y + 6, z + 2);
			}
		}

		for (int y = 8; y < 16; ++y) {
			for (int z = 8; z < 16; ++z) {
				H1_mat.at(y, z) = H(y + 6, z + 6);
			}
		}

		/*
		 ComplexMatrix H1_mat = ComplexMatrix.valueOf(H1_array);

		 H_out = H1_mat.times(DF);
		 */
		//cout << H_1.t() << endl;
		H_out = H1_mat * (DF);

		return H;
}
double** Channel::LTEMIMOCHANNEL(double **Tx1, int tx1length, double **Tx2) {

//////////////////////////////////////////////////////////////////////////////////////put more accurate numbers
	double Medium_type[4][4] = { { 0.84, 0.52, 0.13, 0.08 }, { 0.52, 0.84, 0.08,
			0.13 }, { 0.13, 0.08, 0.84, 0.52 }, { 0.08, 0.13, 0.52, 0.84 } };

	double High_type[4][4] = { { 0.72, 0.45, 0.45, 0.28 }, { 0.45, 0.72, 0.28,
			0.45 }, { 0.45, 0.28, 0.72, 0.45 }, { 0.28, 0.45, 0.45, 0.72 } };

	sqrt_corr_matrix = mat(4, 4, fill::zeros);
	if (corr_type == ("Low")) {
		tx_corr_coeff = 0;
		rx_corr_coeff = 0;
		for (int i = 0; i < 4; ++i) {
			sqrt_corr_matrix(i, i) = 1;
		}
	} else if (corr_type == ("Medium")) {
		tx_corr_coeff = 0.3;
		rx_corr_coeff = 0.9;
		for (int j = 0; j < 4; ++j) {
			for (int i = 0; i < 4; ++i) {
				sqrt_corr_matrix(j, i) = Medium_type[j][i];
			}
		}
	} else {
		tx_corr_coeff = 0.9;
		rx_corr_coeff = 0.9;
		for (int j = 0; j < 4; ++j) {
			for (int i = 0; i < 4; ++i) {
				sqrt_corr_matrix(j, i) = High_type[j][i];
			}
		}
	}
	tx_corr_matrix = mat(2, 2, fill::zeros); //correlation matrix initiallization
	rx_corr_matrix = mat(2, 2, fill::zeros);
	tx_corr_matrix(0, 0) = 1;
	tx_corr_matrix(0, 1) = tx_corr_coeff;
	tx_corr_matrix(1, 0) = tx_corr_coeff;
	tx_corr_matrix(1, 1) = 1;
	rx_corr_matrix(0, 0) = 1;
	rx_corr_matrix(0, 1) = rx_corr_coeff;
	rx_corr_matrix(1, 0) = rx_corr_coeff;
	rx_corr_matrix(1, 1) = 1;
	corr_matrix = kron(tx_corr_matrix, rx_corr_matrix);

	int path_delays[9] = { 0, 30, 70, 90, 110, 190, 410, 0, 0 };
	double path_gains[9] = { 0, -1, -2, -3, -8, -17.2, -20.8, 0, 0 };
	if (channel_type == ("EPA 5Hz")) {

		dopp_freq = 5;
		no_taps = 7;
	} else if (channel_type == ("EVA 5Hz")) {
		path_delays[2] = 150;
		path_delays[3] = 310;
		path_delays[4] = 370;
		path_delays[5] = 710;
		path_delays[6] = 1090;
		path_delays[7] = 1730;
		path_delays[8] = 2510;

		path_gains[1] = -1.5;
		path_gains[2] = -1.4;
		path_gains[3] = -3.6;
		path_gains[4] = -0.6;
		path_gains[5] = -9.1;
		path_gains[6] = -7;
		path_gains[7] = -12;
		path_gains[8] = -16.9;

		dopp_freq = 5;
		no_taps = 9;
	} else if (channel_type == ("EVA 70Hz")) {

		path_delays[2] = 150;
		path_delays[3] = 310;
		path_delays[4] = 370;
		path_delays[5] = 710;
		path_delays[6] = 1090;
		path_delays[7] = 1730;
		path_delays[8] = 2510;

		path_gains[1] = -1.5;
		path_gains[2] = -1.4;
		path_gains[3] = -3.6;
		path_gains[4] = -0.6;
		path_gains[5] = -9.1;
		path_gains[6] = -7;
		path_gains[7] = -12;
		path_gains[8] = -16.9;

		dopp_freq = 70;
		no_taps = 9;
	} else if (channel_type == ("ETU 70Hz")) {
		path_delays[1] = 50;
		path_delays[2] = 120;
		path_delays[3] = 200;
		path_delays[4] = 230;
		path_delays[5] = 500;
		path_delays[6] = 1600;
		path_delays[7] = 2300;
		path_delays[8] = 5000;

		path_gains[0] = -1;
		path_gains[1] = -1;
		path_gains[2] = -1;
		path_gains[3] = 0;
		path_gains[4] = 0;
		path_gains[5] = 0;
		path_gains[6] = -3;
		path_gains[7] = -5;
		path_gains[8] = -7;

		dopp_freq = 70;
		no_taps = 9;
	} else if (channel_type == ("ETU 300Hz")) {
		path_delays[1] = 50;
		path_delays[2] = 120;
		path_delays[3] = 200;
		path_delays[4] = 230;
		path_delays[5] = 500;
		path_delays[6] = 1600;
		path_delays[7] = 2300;
		path_delays[8] = 5000;

		path_gains[0] = -1;
		path_gains[1] = -1;
		path_gains[2] = -1;
		path_gains[3] = 0;
		path_gains[4] = 0;
		path_gains[5] = 0;
		path_gains[6] = -3;
		path_gains[7] = -5;
		path_gains[8] = -7;
		dopp_freq = 300;
		no_taps = 9;
	} else {
		cout << "ERROR : Unknown channel type " << endl;

		dopp_freq = 300;
		no_taps = 9;
	}

///////////////////////////////////////////////////////////////////
	/*
	 cout << tx_corr_matrix.t() << endl;
	 cout << rx_corr_matrix.t() << endl;

	 cout << corr_matrix.t() << endl;
	 cout << sqrt_corr_matrix.t() << endl;
	 cout << (sqrt_corr_matrix * sqrt_corr_matrix).t() << endl;

	 cout << (sqrt_corr_matrix % sqrt_corr_matrix).t() << endl;
	 cout << (sqrt_corr_matrix * sqrt_corr_matrix).t() << endl;*/

///////////////////////////////////////////////////////
	int l = tx1length;

	double f[l * 2];
	for (int k = 0; k < l; ++k) {

		f[k] = fcarry - 59 * 15 * 0.000001 + 15 * 0.000001 * k;
		f[l + k] = fcarry - 59 * 15 * 0.000001 + 15 * 0.000001 * k;

	}

	H = zeros < cx_mat > (2 * l, 2 * l); //zeros(2 * l) + 1i * zeros(2 * l);
	//cout << " Initial H matrix value :" << endl << H.t() << endl;
	//for k = 1:l
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*double tr_1_coeff = 1;						// can be complex...............
	 double tr_2_coeff = 1;*/
	//arma::cx_double test=
	tr_1_coeff = cx_double(1.0, 0.0);
	tr_2_coeff = cx_double(1.0, 0.0);
	tr_1_calc = cx_double(0.0, 0.0);
	tr_2_calc = cx_double(0.0, 0.0);
	A = randn < cx_mat > (4, 1);
	B = sqrt_corr_matrix * A;
	/*
	 cout << " Last A matrix value :" << endl << A.t() << endl;
	 cout << " Last B matrix value :" << endl << B.t() << endl;
	 */

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int k = 0; k < l; ++k) {
		A = randn < cx_mat > (4, 1);	//randn(4, 1) + 1i * randn(4, 1);
		B = sqrt_corr_matrix * A;

		tr_1_coeff = cx_double(1.0, 0.0);
		tr_2_coeff = cx_double(1.0, 0.0);
//cout<<"EXP : "<<tr_1_coeff(0, 0)<<" : "<<exp(tr_1_coeff(0, 0))<<endl;
		// for m =1:no_taps
		for (int m = 0; m < no_taps; ++m) {
			tr_1_calc = cx_double(0.0, 2.0 * PI * f[k] * path_delays[m]);
			tr_2_calc = cx_double(0.0, 2.0 * PI * f[k + 8] * path_delays[m]);

			//	cout << "************** " << m << " " << k <<" EXP : "<<(tr_1_calc(0, 0))<<" : "<<exp(tr_1_calc(0, 0))<< endl;
			tr_1_coeff = tr_1_coeff
					+ sqrt(pow(10.0, (path_gains[m]))) * exp(tr_1_calc);
			tr_2_coeff = tr_2_coeff
					+ sqrt(pow(10.0, (path_gains[m]))) * exp(tr_2_calc);
		}
		//Complex.valueOf(10.0, 0).pow(path_gains[m]).sqrt()
		//.plus(tr_1_calc.exp()).plus(tr_1_coeff);

		H.at(k, k) = tr_1_coeff * B.at(0, 0);        //2 by 2 MIMO --> 4 Paths
		H.at(k, k + l) = tr_2_coeff * B.at(1, 0);
		H.at(k + l, k) = tr_1_coeff * B.at(2, 0);
		H.at(k + l, k + l) = tr_2_coeff * B.at(3, 0);

	}
	//H.print();
	/////////////////

	//////////////////
	/*cout << " tr_1_coeff(0, 0) * B(0, 0) : " << tr_1_coeff(0, 0) << " * "
	 << B(0, 0) << " = " << tr_1_coeff(0, 0) * B(0, 0) << endl;
	 cout << " Last A matrix value :" << endl << A.t() << endl;
	 cout << " Last B matrix value :" << endl << B.t() << endl;

	 cout << H.t() << endl;*/
	///////////////////////////////////////////////////////////////
	noise = randn < cx_mat > (2 * l, 1);

	/*double **Tx_test;
	Tx_test = Tx1;
	//FFT_12 fft = FFT_12();
	Tx_test = fft.getfftw(12, Tx_test, 12);
	Tx_test = fft.getfftw(12, Tx1, 12);*/
//	Tx1 = fft.getfftw(12, Tx1, 12);
//	Tx2 = fft.getfftw(12, Tx2, 12);


	for (int var = 0; var < 12; ++var) {
		Tx.at(var, 0) = cx_double(Tx1[0][var], Tx1[1][var]);
		Tx.at(var + 12, 0) = cx_double(Tx2[0][var], Tx2[1][var]);
		//cout<<Tx1[0][var]<<"+"<< Tx1[1][var]<<"i * ";
	}
	Tx.print("This is channel tx: ");
	//cout<<"This is TX vector "<< Tx.t()<<endl;

///////////////////////////////////////////////////////////////////////////////////////////
	R = H * Tx;	// + sigma * noise;

	//cout<<"This is R vector "<< R.t()<<endl;

	for (int var = 0; var < 24; ++var) {
		Rx[0][var] = R.at(var, 0).real();
		Rx[1][var] = R.at(var, 0).imag();

		//cout<<Rx[0][var]<<"+"<<Rx[1][var]<<"i \t";
	}
//	cout << "this is RX b: " << Rx[0][0] << endl;
//	Rx = ifft.getifftw_channel(12, Rx, 24);

//	cout << "this is RX a: " << Rx[0][0] << endl;
	//Rx = ifft(Rx);
//	arma::cx_mat Rx1 = ifft(R(1:8));
//	arma::cx_mat Rx2 = ifft(R(9:16));
	///////////////////////////////////////////////////////////////////////////////////////

	//arma::cx_mat H_1 = zeros<cx_mat>(16, 16);

	//for y=1:8,
	for (int y = 0; y < 8; ++y) {
		for (int z = 0; z < 8; ++z) {
			//		H_1(y, z) = H(y + 2, z + 2);
			H1_mat.at(y, z) = H(y + 2, z + 2);
		}
	}

	for (int y = 0; y < 8; ++y) {
		for (int z = 8; z < 16; ++z) {
			H1_mat.at(y, z) = H(y + 2, z + 6);
		}
	}

	for (int y = 8; y < 16; ++y) {
		for (int z = 0; z < 8; ++z) {
			H1_mat.at(y, z) = H(y + 6, z + 2);
		}
	}

	for (int y = 8; y < 16; ++y) {
		for (int z = 8; z < 16; ++z) {
			H1_mat.at(y, z) = H(y + 6, z + 6);
		}
	}

	/*
	 ComplexMatrix H1_mat = ComplexMatrix.valueOf(H1_array);

	 H_out = H1_mat.times(DF);
	 */
	//cout << H_1.t() << endl;
	H_out = H1_mat * (DF);

	return Rx;
}
