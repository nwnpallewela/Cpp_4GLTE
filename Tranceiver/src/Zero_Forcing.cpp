//============================================================================
// Name        : Zero_Forcing.cpp
// Author      : Nuwan Pallewela
// Version     :
// Copyright   : Feel free
// Description :
//============================================================================

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

#include"demapper.h"
#include"modulator.h"
#include"fft.h"
#include"fft12.h"
#include"ifft.h"
#include"dftarray.h"
#include"channel.h"
#include"demodulator.h"
#include"equalizer.h"
#include"decoder.h"
#include"ifft_zeroforcing.h"
#include"ML_zeroforcing.h"

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

/////////////////////////////////////////////////////////////////////////////////////
void dataGen(int size);
string readData(int size);
string turboencode(string data, int blocksize);	// data,blocksize of data
string interleave(string data_in, int K);
string multiplex(string code1, string code2);
string encode(string DATA, int K);
int error_calc(string rawdata_in, int* decoded_data, int count);

//////////////////////////////////////////////////////////////////////////////////////
int main() {
	for (int var = 0; var < 1; ++var) {

		srand((unsigned) time(0));
		cout << endl << endl
				<< "***********************************************************************"
				<< endl
				<< "***********************************************************************"
				<< endl
				<< "   This is a test of LTEA C++ Transceiver implementation Zero Forcing "
				<< endl
				<< "***********************************************************************"
				<< endl
				<< "***********************************************************************"
				<< endl << endl;

		int datafilesize = 1;
		int data_block_size = 40;
		int decoder_num = 2;
		int error = 0;
		////////////////////////////////////////////////////////////////////////////////// Enter data fle size in console
		/*	while (datafilesize % 40 != 0) {
		 cout << endl << "Enter data file size ( X 320 ): ";

		 cin >> datafilesize;

		 cout << endl;
		 if (datafilesize % 320 != 0) {
		 cout << endl << "ERROR : The size should be a multiple of 320 "
		 << endl;
		 }

		 }*/
		datafilesize = 320 * 10000;
		//////////////////////////////////////////////////////////////////////////////////
		int blocksize = 40;
		int starttime = 0;
		int endtime = 0;
		int len = 0;
		//timeval t1, t2;

		//generate data file
		double data_write_read_time = 0;
		double encode_time_full = 0;
		double modulate_time_full = 0;
		double fft_time_full = 0;
		double ifft_time_full = 0;
		double multiplex_time_full = 0;
		double decoder_time_full = 0;
		//double decoder_time_full = 0;
		double channel_time_full = 0;
		double demapper_time_full = 0;
		double equalizer_time_full = 0;
		double demodulator_time_full = 0;
		double getH_time_full = 0.0;

		double total_time_without_dataread = clock();
		string rawdata = "";
		string encodeddata = "";
		string received_data = "";
		int antenna1_count = 0;
		int antenna2_count = 0;
		double **Rx = new double*[2];

		for (int i = 0; i < 2; i++) {
			Rx[i] = new double[24];

			for (int j = 0; j < 24; j++) {
				Rx[i][j] = 0;
			}
		}

		double **y_out = new double*[2];
		for (int i = 0; i < 2; i++) {
			y_out[i] = new double[16];

			for (int j = 0; j < 16; j++) {
				y_out[i][j] = 0;
			}
		}

		double **antenna1 = new double*[2];
		double **antenna2 = new double*[2];
		for (int i = 0; i < 2; i++) {
			antenna1[i] = new double[132];
			antenna2[i] = new double[132];
			for (int j = 0; j < (132); j++) {
				antenna1[i][j] = 0;
				antenna2[i][j] = 0;
			}
		}
		double **antenna1_ch = new double*[2];
		double **antenna2_ch = new double*[2];
		for (int i = 0; i < 2; i++) {
			antenna1_ch[i] = new double[12];
			antenna2_ch[i] = new double[12];
			for (int j = 0; j < 12; j++) {
				antenna1_ch[i][j] = 0;
				antenna2_ch[i][j] = 0;
			}
		}
		double** Tx1;
		double** Tx2;

		double **Tx_1;
		double **Tx_2;

		double **fft_data;
		double **modulateddata;
		double **ifftdata = new double*[2];
		for (int i = 0; i < 2; i++) {
			ifftdata[i] = new double[(264)];
			for (int j = 0; j < (264); j++) {
				ifftdata[i][j] = 0;
			}
		}
		double *y_ = new double[32];
		FFT fft_ob = FFT(176);
		IFFT ifft_ob = IFFT(176);
		Modulator modu = Modulator(176);	//encodeddata.length() / 6;

		FFT_12 fft;
		IFFT_CH ifft;
		cx_mat H;
		cx_mat noise;
		cx_mat Tx = zeros<cx_mat>(24, 1);

		cx_mat R;
		DeMapper DM = DeMapper();
		//Equalizer Eq = Equalizer();

		IFFT_Z IFFT_z = IFFT_Z(16);
		ML_Z ML = ML_Z();

		Demodulator Dmod = Demodulator();
		arma::mat Y = mat(32, 1, fill::zeros);
		arma::mat S = mat(32, 1, fill::zeros);
		Decoder dec = Decoder();
		Channel ch = Channel("High", "EPA 5Hz", 5, dec.getSigma());
		cout << dec.getSigma() << endl;
		double lc = 2.5;
		int iterations = 6;
		double *LLR = new double[40];	//=  double[40];
		double *LLR1 = new double[40];
		double *LLR2 = new double[40];
		/* double A1[][]=new double[17][44];
		 double A2[][]=new double[17][44];*/
		//	int I1[][]=new int[40][2];
		int **I2 = new int*[2];
		I2[0] = new int[40];
		I2[1] = new int[40];
		double *leuk1 = new double[40];
		double *leuk2 = new double[40];
		double *luk = new double[40];
		double *R1 = new double[40];
		double *R2 = new double[40];
		int *y1 = new int[40];
		int *y2 = new int[40];
		double *LEUK = new double[40];
		int* decode;
		decode = new int[40];
		int* decoded_data;
		decoded_data = new int[320];
		for (int i_size = 0; i_size < datafilesize; i_size = i_size + 320) {

			starttime = clock();
			dataGen(320);
			rawdata = readData(320);				//read from file
			endtime = clock();
			data_write_read_time = data_write_read_time + (endtime - starttime);

			/////////////////////////////////////////////////////////////////////////////////////encode data using turbo encoder

			starttime = clock();
			encodeddata = turboencode(rawdata, blocksize);
			endtime = clock();
			encode_time_full = (endtime - starttime) + encode_time_full;
			////////////////////////////////////////////////////////////////////////////////////print encoded data, time and length

			///////////////////////////////////////////////////////////////////////////////////modulating the data

			starttime = clock();
			modulateddata = modu.modulate64QAM(encodeddata);
			endtime = clock();

			modulate_time_full = (endtime - starttime) + modulate_time_full;
			//////////////////////////////////////////////////////////////////////////////////fft
			starttime = clock();
			//double **fft_data;
			fft_data = modulateddata;

			///////////////////////////////////////////////////////////////////////////////////////fft mapping and dft
			starttime = clock();
			len = (int) (encodeddata.length() / 6);
			//FFT fft_ob=  FFT(len);

			fft_data = fft_ob.getfftw(8, modulateddata, len);
			endtime = clock();

			fft_time_full = (endtime - starttime) + fft_time_full;
			//////////////////////////////////////////////////////////////////////////////////
			starttime = clock();
			ifftdata = ifft_ob.getifftw(12, fft_data,
					(int) (encodeddata.length() / 6));
			endtime = clock();

			ifft_time_full = (endtime - starttime) + ifft_time_full;
			///////////////////////////////////////////////////////////////////////////////////////divide data to two antennas

			antenna1_count = 0;
			antenna2_count = 0;
			starttime = clock();
			for (int j = 0; j < (encodeddata.length() / 4); j++) {
				if ((j % 24) < 12) {
					antenna1[0][antenna1_count] = ifftdata[0][j];
					antenna1[1][antenna1_count] = ifftdata[1][j];
					antenna1_count++;

				} else {
					antenna2[0][antenna2_count] = ifftdata[0][j];
					antenna2[1][antenna2_count] = ifftdata[1][j];
					antenna2_count++;
				}
			}
			endtime = clock();

			multiplex_time_full = (endtime - starttime) + multiplex_time_full;

			///////////////////////////////////////////////////////////////////////////////// channel execution
			// ////////////////////////////////////////////////////////////////////////////////////////////////////
			//
			//
			// Each iteration of the transmitter send 11 antenna outs from each
			// antenna
			// in sets of 12 complex numbers. So channel and receiver have to
			// run 11 iterations
			// to process all transmitted data
			//
			// ////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < 2; i++) {

				for (int j = 0; j < 12; j++) {
					antenna1_ch[i][j] = 0;
					antenna2_ch[i][j] = 0;
				}
			}

			received_data = "";

			for (int i = 0; i < 11; i++) {
				// count1 = 0;

				starttime = clock();
				for (int j = 0; j < 12; j++) {
					antenna1_ch[0][j] = antenna1[0][12 * i + j];
					antenna1_ch[1][j] = antenna1[1][12 * i + j];

					antenna2_ch[0][j] = antenna2[0][12 * i + j];
					antenna2_ch[1][j] = antenna2[1][12 * i + j];
					// count1++;
					//cout<<antenna1_ch[0][j]<<" + "<<antenna1_ch[1][j]<<"i          ---------          "<<antenna2_ch[0][j]<<" + "<<antenna2_ch[1][j]<<endl;
				}

				//	Rx = ch.LTEMIMOCHANNEL(antenna1_ch, 12, antenna2_ch);
				/////////////////////////////////////////////////////////////////////////////////////////
				H = ch.LTEMIMOCHANNEL_mat();

				noise = randn<cx_mat>(2 * 12, 1);
//noise.print();
				Tx_2 = fft.getfftw(12, antenna2_ch, 12);

				Tx_1 = fft.getfftw(12, antenna1_ch, 12);

				for (int var = 0; var < 12; ++var) {
					Tx.at(var, 0) = cx_double(Tx_1[0][var], Tx_1[1][var]);

				}
				Tx_1 = fft.getfftw(12, antenna1_ch, 12);

				Tx_2 = fft.getfftw(12, antenna2_ch, 12);

				for (int var = 0; var < 12; ++var) {
					Tx.at(var + 12, 0) = cx_double(Tx_2[0][var], Tx_2[1][var]);
				}
//Tx.print();
				//	H.print("This is h");
				///////////////////////////////////////////////////////////////////////////////////////////
				R = (H * Tx);
				//	R.print("This is R before");
				noise = dec.getSigma() * noise;
				R = R + noise;
				/*noise.print("this is noise");*/
//noise=(dec.getSigma() * noise);
//noise.print("thisis noise*sigma");
			//	R.print("This is R");

				for (int var = 0; var < 24; ++var) {
					Rx[0][var] = R.at(var, 0).real();
					Rx[1][var] = R.at(var, 0).imag();

				}
				Rx = ifft.getifftw_channel(12, Rx, 24);
				endtime = clock();
				channel_time_full = channel_time_full + (endtime - starttime);

				//////////////////////////////////////////////////////////////////////////////////////////
//arma::inv(H);
				starttime = clock();
				//H.print("this is H");
				//(H*arma::inv(H)).print("H*Hinv");
				y_ = DM.get_demapped_rx_ZF(Rx, arma::inv(H));
				endtime = clock();
				demapper_time_full = demapper_time_full + (endtime - starttime);
				//	cout<<"******************"<<endl;
				starttime = clock();

			/*	cout << "this is Y" << endl;
				for (int var = 0; var < 32; ++var) {
					cout << y_[var] << endl;
				}*/
				//Eq.genH(ch.getHout());
				endtime = clock();
				getH_time_full = getH_time_full + (endtime - starttime);
				starttime = clock();
				for (int var = 0; var < 32; ++var) {
					Y.at(var, 0) = y_[var];
				}

				y_out = IFFT_z.getifftw(8, Y, 16);
				//	cout << "this is Ydash" << endl;
				//	starttime = clock();
				/* for (int var = 0; var < 16; ++var) {
				 cout << y_out[0][var] <<" + "<<y_out[1][var] <<"i"<< endl;
				 }*/
				y_ = ML.Decision(y_out);
				/*cout << "this is Y" << endl;

				for (int var = 0; var < 32; ++var) {
				 cout << y_[var] << endl;
				 }*/
				//y_ = Eq.GetLSD_Y(Y);
				endtime = clock();
				equalizer_time_full = equalizer_time_full
						+ (endtime - starttime);

				starttime = clock();

				received_data = received_data + Dmod.runDemodulator(y_);
				endtime = clock();
				demodulator_time_full = demodulator_time_full
						+ (endtime - starttime);
			}

			int count = 0;

			starttime = clock();

			for (int i1 = 0; i1 < encodeddata.length(); i1 = i1 + 132) {

				dec.decoder_log_map_it(received_data.substr(i1, i1 + 132), 1,
						data_block_size, luk);

				LLR1=dec.getLLR1();
			/*	cout<<"This is LLR1"<<endl;
				for (int var = 0; var < 40; ++var) {
					cout<<LLR1[var]<<endl;
				}*/



				leuk1 = dec.get_leuk();
				y1 = dec.getY();
				dec.decoder_log_map_it(received_data.substr(i1, i1 + 132), 2,
						data_block_size, leuk1);
				/*LLR2=dec.getLLR1();
				cout<<"This is LLR2"<<endl;
								for (int var = 0; var < 40; ++var) {
									cout<<LLR2[var]<<endl;
								}*/
				I2 = dec.get_interleave_table();
				leuk2 = dec.get_leuk();
				R2 = dec.getR();
				y2 = dec.getY();
				LEUK = leuk2; //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				lc=dec.getLC();
				for (int n = 3; n <= iterations; ++n) { //   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

					if ((n % 2) == 1) {
						//  for i=1:1:data_block_size
						for (int i = 0; i < data_block_size; ++i) {
							LLR[i] = R1[i] + LEUK[i];
							LEUK[i] = LLR[i] - LEUK[i] - 2 * lc * y1[2 * i];
						}
					} else {
						for (int i = 0; i < data_block_size; ++i) {
							LLR[I2[i][1]] = R2[i] + LEUK[I2[i][1]]; //%%%%%%%%%%%%%%%%%%%%%%%
							LEUK[I2[i][1]] = LLR[I2[i][1]] - LEUK[I2[i][1]]
									- 2 * lc * y2[2 * i]; //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						}

					}

				} //                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				for (int i = 0; i < 40; i++)
				{
				//	cout<<"*  "<<LLR[i]<<endl;
					if (LLR[i] > 0) {
						decode[i] = 1;
					} else {
						decode[i] = 0;
					}

				}

				for (int j = 0; j < 40; ++j) {
					decoded_data[count] = decode[j];
					count++;
				}

			}
			endtime = clock();

			decoder_time_full = decoder_time_full + (endtime - starttime);
//cout<<rawdata<<endl<<decoded_data<<;
			error = error + error_calc(rawdata, decoded_data, rawdata.length());

		}
		cout << endl << "**************************************** " << endl;

		total_time_without_dataread = clock() - total_time_without_dataread
				- data_write_read_time;

		cout << endl << "Data file size           = " << datafilesize << endl;
		cout << "data read write time     = " << data_write_read_time << " us"
				<< endl << endl;
		cout << endl << "turbo encode time        = " << encode_time_full
				<< " us" << endl;
		cout << "modulate time            = " << modulate_time_full << " us"
				<< endl;
		cout << "fft time                 = " << fft_time_full << " us" << endl;
		cout << "ifft time                = " << ifft_time_full << " us"
				<< endl;
		cout << "Spatial multiplexing time= " << multiplex_time_full << " us"
				<< endl;

		cout << "Channel time             = " << channel_time_full << " us"
				<< endl;
		cout << "Get H time               = " << getH_time_full << " us"
				<< endl;
		cout << "Demapper time            = " << demapper_time_full << " us"
				<< endl;
		cout << "Equalizer time           = " << equalizer_time_full << " us"
				<< endl;
		cout << "Demodulator time         = " << demodulator_time_full << " us"
				<< endl;
		cout << "Decoder time             = " << decoder_time_full << " us"
				<< endl;

		cout << endl << "total time               = "
				<< (encode_time_full + modulate_time_full + fft_time_full
						+ ifft_time_full + multiplex_time_full
						+ demapper_time_full + equalizer_time_full
						+ demodulator_time_full + decoder_time_full) / 1000
				<< " ms" << endl;

		cout << "total_time_without_dataread = "
				<< total_time_without_dataread / 1000 << " ms" << endl;

		cout << "Error count           = " << error << " / " << datafilesize
				<< endl;

		cout << endl << "finished transmitter " << endl;

////////////////////////////////////////////////////////////////////////
	}
}

string turboencode(string data, int blocksize) {
////////////////////////////////////////////////////////////////variables decleration
	string encodeddata = encode(data, blocksize);

///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
	return encodeddata;
}

string encode(string DATA, int K) {

	string code1 = "";
	string code2 = "";
	string code_word = "";
	string data_atr = "";
	string interleaved_data = "";
// /////////////////////////////////////////////////////////

	int j = 0;

	int D0 = 0;
	int D1 = 0;
	int D2 = 0;
	int D3 = 0;

	int tD0 = 0;
	int tD1 = 0;
	int tD2 = 0;
	int tD3 = 0;
	int j1 = 0;

	int D02 = 0;
	int D12 = 0;
	int D22 = 0;
	int D32 = 0;

	int tD02 = 0;
	int tD12 = 0;
	int tD22 = 0;
	int tD32 = 0;
	int j12 = 0;

// System.out.println(data.length()-K);
	for (j = 0; j <= (DATA.length() - K); j = j + K) {
		int i1 = 0;
		int i12 = 0;

		data_atr = DATA.substr(j, K);

//data_atr = data_atr.concat("2");
		for (int i = 0; i < K; i++) {
			tD0 = D0;
			tD1 = D1;
			tD2 = D2;
			tD3 = D3;
			D0 = data_atr[i] - 48;
			D3 = tD2;
			D2 = tD1;
			D1 = ((D0 + tD2 + tD3) % 2);

			//cout << D0 << " " << D1 << " " << D2 << " " << (D3) << endl;
			if (D0 == 0) {
				code1 = code1 + "0";
			} else {
				code1 = code1 + "1";
			}
			if (((D1 + tD1 + tD3) % 2) == 0) {
				code1 = code1 + "0";
			} else {
				code1 = code1 + "1";
			}

			i1 = 2 * i;
		}

		for (int k = 0; k < 3; k++) {
			tD0 = D0;
			tD1 = D1;
			tD2 = D2;
			tD3 = D3;
			D0 = ((tD2 + tD3) % 2);
			D3 = tD2;
			D2 = tD1;
			D1 = ((D0 + tD2 + tD3) % 2);

			//D1 = tD0;

			if (D0 == 0) {
				code1 = code1 + "0";
			} else {
				code1 = code1 + "1";
			}
			if (((D1 + tD1 + tD3) % 2) == 0) {
				code1 = code1 + "0";
			} else {
				code1 = code1 + "1";
			}

		}

		interleaved_data = interleave(data_atr, K);

		for (int i = 0; i < K; i++) {

			tD02 = D02;
			tD12 = D12;
			tD22 = D22;
			tD32 = D32;
			D02 = interleaved_data[i] - 48;
			D32 = tD22;
			D22 = tD12;
			D12 = ((D02 + tD22 + tD32) % 2);
			if (D02 == 0) {
				code2 = code2 + "0";
			} else {
				code2 = code2 + "1";
			}
			if (((D12 + tD12 + tD32) % 2) == 0) {
				code2 = code2 + "0";
			} else {
				code2 = code2 + "1";
			}

			i12 = 2 * i;
		}
		for (int k = 0; k < 3; k++) {
			tD02 = D02;
			tD12 = D12;
			tD22 = D22;
			tD32 = D32;
			D02 = ((tD22 + tD32) % 2);
			D32 = tD22;
			D22 = tD12;

			D12 = ((D02 + tD22 + tD32) % 2);
			//	D12 = (tD02);

			if (D02 == 0) {
				code2 = code2 + "0";
			} else {
				code2 = code2 + "1";
			}
			if (((D12 + tD12 + tD32) % 2) == 0) {
				code2 = code2 + "0";
			} else {
				code2 = code2 + "1";
			}

		}

		//	cout << "code1 " << code1 << endl;
		//	cout << "code2 " << code2 << endl;
		code_word = code_word
				+ (multiplex(code1.substr(j1, i1 + 2),
						code2.substr(j12, i12 + 2)));
// System.out.println(code_word.length());
		code_word = code_word + (code1.substr(j1 + i1 + 2, 6));
// System.out.println(code_word.length());
		code_word = code_word + (code2.substr(j12 + i12 + 2, 6));

		j1 = j1 + 2 * K + 6;
		j12 = j12 + 2 * K + 6;

	}

// //////////////////////////////////////////////////////////
// System.out.println(code_word.length());
	return code_word;
}

string multiplex(string code1, string code2) {

// System.out.println(code1.length());
	string arr = "";

	for (int i = 0; i < code1.length(); i = i + 2) {

		arr = arr + (code1[i]);

		arr = arr + (code1[i + 1]);
		arr = arr + (code2[i + 1]);

	}

	return arr;
}

string interleave(string data_in, int K) {
	int f1 = 3;
	int f2 = 10;
	int index;

	string data_out = "";

	for (int i = 0; i < K; i++) {
		index = (((f1 * (i) + f2 * (i) * (i))) % (K));

		data_out = data_out + data_in[index];
		//cout << " interleave " << i << " " << index << endl;
	}
//cout << endl;
	return data_out;
}

void dataGen(int size) {

	ofstream out("data");
	int i = 0;
	if (!out) {
		cout << "Cannot open data file.\n";

	}

	for (int i = 0; i < size; i++) {
		//cout<<rand()<<endl;
		if ((rand() % 2) == 0) {
			out << "0";
		} else {
			out << "1";
		}

	}

	out.close();

}

string readData(int size) {
	ifstream in("data"); // input
	//ifstream in("ones");
	if (!in) {
		cout << "Cannot open data file.\n";
		return "Woo";
	}
	char item[size];

	in >> item;

	in.close();
	return item;
}

int error_calc(string rawdata_in, int* decoded_data, int count) {
	int error = 0;

	if (rawdata_in.length() == count) {
		for (int i = 0; i < rawdata_in.length(); ++i) {
		//	cout<<decoded_data[i]<<" -- "<<((int) rawdata_in[i] - 48)<<endl;
			if (((int) rawdata_in[i] - 48) != decoded_data[i]) {
				error++;
			}
		}
		return error;
	} else {
		// cout << " ERROR : Bit counts are different : " <<
		// rawdata_in.length()
		// << " " << count << endl;
		return 677;
	}
}
