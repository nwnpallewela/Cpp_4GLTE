//============================================================================
// Name        : 4GLTECpp.cpp
// Author      : Nuwan Pallewela
// Version     : transmitter fftw3 //
// Copyright   : Feel free
// Description : Hello World in C++, Ansi-style
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

using namespace std;

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

/////////////////////////////////////////////////////////////////////////////////////
void dataGen(int size);
string readData(int size);
string turboencode(string data, int blocksize);	// data,blocksize of data
string interleave(string data_in, int K);
string multiplex(string code1, string code2);
string encode(string DATA, int K);
int findindex(string code);
double **modulate64QAM(string encodeddata);
double **getfftw(int N, double **modulateddata, int len);
double **getifftw(int N, double **modulateddata, int len);

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

int main() {
	srand((unsigned) time(0));
	cout << endl << endl << "*************************************************"
			<< endl << "*************************************************"
			<< endl << "This is a test of LTEA C++ Transmitter implementation "
			<< endl << "*************************************************"
			<< endl << "*************************************************"
			<< endl << endl;

	int datafilesize = 1;
	////////////////////////////////////////////////////////////////////////////////// Enter data fle size in console
	while (datafilesize % 40 != 0) {
		cout << endl << "Enter data file size ( X 40 ): ";

		cin >> datafilesize;

		cout << endl;
		if (datafilesize % 40 != 0) {
			cout << endl << "ERROR : The size should be a multiple of 40 "
					<< endl;
		}

	}
	//////////////////////////////////////////////////////////////////////////////////
	int blocksize = 40;
	int starttime = 0;
	int endtime = 0;
	//timeval t1, t2;

	starttime = clock();
	dataGen(datafilesize);									//generate data file
	string rawdata = readData(datafilesize);				  	//read from file
	endtime = clock();

	//cout << "input data stream : " << rawdata << endl;

	cout << "data read write time     = " << (endtime - starttime) << " us"
			<< endl << endl;
	/////////////////////////////////////////////////////////////////////////////////////encode data using turbo encoder

	starttime = clock();

	string encodeddata = turboencode(rawdata, blocksize);
	endtime = clock();

	////////////////////////////////////////////////////////////////////////////////////print encoded data, time and length

	//cout << "Encoded data stream : " << endl << endl << encodeddata << endl;

	cout << endl << "Encoded data length : " << encodeddata.length() << endl;
	//cout <<endl<< "turbo encode time        = " << (endtime - starttime) << " us" << endl;
	int encode_time = (endtime - starttime);
	///////////////////////////////////////////////////////////////////////////////////modulating the data
	starttime = clock();
	double **modulateddata = modulate64QAM(encodeddata);
	endtime = clock();
//	cout << "modulation time          = " << (endtime - starttime) << " us" << endl;
	int modulate_time = (endtime - starttime);
	//////////////////////////////////////////////////////////////////////////////////fft
	starttime = clock();
	double **fft_data;
	fft_data = modulateddata;
	double **ifftdata = new double*[2];
	for (int i = 0; i < 2; i++) {
		ifftdata[i] = new double[(encodeddata.length() / 4)];
		for (int j = 0; j < (encodeddata.length() / 4); j++) {
			ifftdata[i][j] = 0;
		}
	}
///////////////////////////////////////////////////////////////////////////////////////fft mapping and dft
	starttime = clock();
	fft_data = getfftw(8, modulateddata, (int) (encodeddata.length() / 6));
	endtime = clock();
	//cout << "fft time                 = " << (endtime - starttime) << " us" << endl;
	int fft_time = (endtime - starttime);
	//////////////////////////////////////////////////////////////////////////////////
	starttime = clock();
	ifftdata = getifftw(12, fft_data, (int) (encodeddata.length() / 6));
	endtime = clock();
	//cout << "ifft time                = " << (endtime - starttime) << " us" << endl;
	int ifft_time = (endtime - starttime);
///////////////////////////////////////////////////////////////////////////////////////divide data to two antennas

	double **antenna1 = new double*[2];
	double **antenna2 = new double*[2];
	for (int i = 0; i < 2; i++) {
		antenna1[i] = new double[(encodeddata.length() / 8)];
		antenna2[i] = new double[(encodeddata.length() / 8)];
		for (int j = 0; j < (encodeddata.length() / 8); j++) {
			antenna1[i][j] = 0;
			antenna2[i][j] = 0;
		}
	}

	int antenna1_count = 0;
	int antenna2_count = 0;
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
	//cout << "spacial multiplexing time= " << (endtime - starttime) << " us"
	//		<< endl;
	int multiplex_time = (endtime - starttime);

	//cout <<endl<< "total time               = "
	//		<< encode_time + modulate_time + fft_time + ifft_time
	//				+ multiplex_time << " us" << endl;

	//////////////////////////////////////////////////////////////////////////////////////printing fft data

	/*	for (int i = 0; i < 2; i++) {

	 for (int j = 0; j < (encodeddata.length() / 6); j++) {
	 cout << fft_data[i][j] << " ";
	 }
	 cout << endl;
	 }*/
/////////////////////////////////////////////////////////////////////////////////////////printing ifft data
//	for (int i = 0; i < 2; i++) {
//		for (int j = 0; j < (encodeddata.length() / 4); j++) {
//			if (ifftdata[i][j] > 0.0005 || ifftdata[i][j] < -0.0005)
//				cout << ifftdata[i][j] << " ";
//			else
//				cout << "0" << " ";
//		}
//		cout << endl;
//	}
	/////////////////////////////////////////////////////////////////////////////////////printing antenna data
	/*	cout << endl << "    : antenna 1 data :            : antenna  2  data :  "
	 << endl << endl;

	 for (int j = 0; j < (encodeddata.length() / 8); j++) {

	 printf("%9.5f , %9.5f i   :    %9.5f , %9.5f i \n", antenna1[0][j],
	 antenna1[1][j], antenna2[0][j], antenna2[1][j]);

	 }*/

	cout << endl << "Data file size           = " << datafilesize << endl;
	cout << endl << "turbo encode time        = " << encode_time << " us"
			<< endl;
	cout << "modulate time            = " << modulate_time << " us" << endl;
	cout << "fft time                 = " << fft_time << " us" << endl;
	cout << "ifft time                = " << ifft_time << " us" << endl;
	cout << "spacial multiplexing time= " << (endtime - starttime) << " us"
			<< endl;
	cout << endl << "total time               = "
			<< encode_time + modulate_time + fft_time + ifft_time
					+ multiplex_time << " us" << endl;
	cout << endl << "finished transmitter " << endl;
////////////////////////////////////////////////////////////////////////

}

double **getfftw(int N, double **modulateddata, int len) {

	fftw_complex in[N], out[N]; /* double [2] */
	fftw_plan p;

	double **fftdata = new double*[2];
	for (int i = 0; i < 2; i++) {
		fftdata[i] = new double[len];
		for (int j = 0; j < len; j++) {
			fftdata[i][j] = 0;
		}
	}

	for (int j = 0; j < len; j = j + 8) {
		for (int i = 0; i < N; i++) {
			in[i][0] = modulateddata[0][j + i];
			in[i][1] = modulateddata[1][j + i];
		}

		/* forward Fourier transform, save the result in 'out' */
		p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		for (int i = 0; i < N; i++) {
			fftdata[0][i + j] = out[i][0];
			fftdata[1][i + j] = out[i][1];
		}

	}
	fftw_destroy_plan(p);

	return fftdata;

}

double **getifftw(int N, double **fftdata, int len) {
	fftw_complex in2[N], out[N]; /* double [2] */
	fftw_plan q;
	double **ifftdata = new double*[2];
	for (int i = 0; i < 2; i++) {
		ifftdata[i] = new double[len * 3 / 2];
		for (int j = 0; j < len; j++) {
			ifftdata[i][j] = 0;
		}
	}
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

		q = fftw_plan_dft_1d(N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(q);

		/* normalize */
		for (int i = 0; i < N; i++) {
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
	fftw_destroy_plan(q);
	fftw_cleanup();

	return ifftdata;
}

string turboencode(string data, int blocksize) {
////////////////////////////////////////////////////////////////variables decleration
	string encodeddata = encode(data, blocksize);

///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
	return encodeddata;
}

double **modulate64QAM(string encodeddata) {

	const int real[64] = { -7, -7, -7, -7, -7, -7, -7, -7, -5, -5, -5, -5, -5,
			-5, -5, -5, -1, -1, -1, -1, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3,
			-3, -3, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1,
			1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3 };
	const int complex[64] = { -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5,
			1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7,
			-5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3,
			7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3 };
//////////////////////////////////////////////////////////////////////////////////////////////////////////

	double **modulateddata = new double*[2];
	for (int i = 0; i < 2; i++) {
		modulateddata[i] = new double[encodeddata.length() / 6];
		for (int j = 0; j < encodeddata.length() / 6; j++) {
			modulateddata[i][j] = 0;
		}
	}
	int index;
	int count = 0;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j <= (encodeddata.length() - 6); j = j + 6) {

		index = findindex(encodeddata.substr(j, 6));

		modulateddata[0][count] = real[index];
		modulateddata[1][count] = complex[index];
		count++;
//	cout<<encodeddata.substr(j, 6)<<" "<<real[index]<<" "<<complex[index]<<endl;
	}

	return modulateddata;
}

int findindex(string code) {

	return ((code[0] - 48) + (code[2] - 48) * 2 + (code[2] - 48) * 4
			+ (code[3] - 48) * 8 + (code[4] - 48) * 16 + (code[5] - 48) * 32);
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
			out << "1";
		} else {
			out << "0";
		}

	}

	out.close();

}

string readData(int size) {
	ifstream in("data"); // input
	if (!in) {
		cout << "Cannot open data file.\n";
		return "Woo";
	}
	char item[size];

	in >> item;

	in.close();
	return item;
}
