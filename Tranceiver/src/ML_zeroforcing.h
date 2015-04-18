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


class ML_Z {

private:
	double *y_;

public:
	 ML_Z (){
		y_=new double[32];

	}
	 double* Decision(double** yout){

		for (int i = 0; i < 16; i++) {
			y_[i]=(yout[0][i]+1.0)/2.0;
			y_[i+16]=(yout[1][i]+1.0)/2.0;
		}
		for (int i = 0; i < 32; i++){
			y_[i]=round(y_[i]);
			if(y_[i]>4){
				y_[i]=4;
			}else if(y_[i]<-3){
				y_[i]=-3;
			}
			y_[i]=2*y_[i]-1;
		}

		/*for (int i = 0; i < 16; i++) {
				cout<<	y_[i]<<"  --  "<<(yout[0][i])<<endl;
					cout<<y_[i+16]<<"  --  "<<(yout[1][i])<<endl;
				}*/

	/*	for (int i = 0; i < 16; i++) {
			System.out.println(yout[0][i]+" + "+yout[1][i]+" ---- "+y_[i]+" + "+y_[i+16]+"i");
		}*/

		return y_;
	}

};
