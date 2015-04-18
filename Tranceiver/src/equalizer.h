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

#include"lsdtree.h"

using namespace std;
using namespace arma;

class Equalizer {
private:
	mat H;
	mat Q;
	mat R;
	mat Rs;
	Node minNode;
	LSDTree Tree;


	double *y;

public:
	void genH(cx_mat H1);
	void RS(mat S);
	double* GetLSD_Y(mat Y);
	mat getRS(mat S);
	mat getQtRS(mat S);
	Equalizer() {
		Tree = LSDTree(1000000.0, 32); //R and K
		H = mat(32, 32, fill::zeros);
		Q = mat(32, 32, fill::zeros);
		R = mat(32, 32, fill::zeros);
		Rs = mat(32, 1, fill::zeros);
		minNode=Node();
		y=new double[32];
	}

};
void Equalizer::genH(cx_mat H1) {

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {

			H.at(i, j) = H1.at(i, j).real();
			H.at(i + 16, j + 16) = H1.at(i, j).real();
			H.at(i + 16, j) = H1.at(i, j).imag();
			H.at(i, j + 16) = (-1) * H1.at(i, j).imag();
		}
	}

	//H1.print();

}

double* Equalizer::GetLSD_Y(mat Y) {

	//QRDecomposition QRD = new QRDecomposition(H);
	/*

	arma::qr(Q, R, H);
//	(H - Q * R).print("This is H-QR");
	Q=QRD.getQ();
	 R=QRD.getR();

	System.out.println("Y row: ");
	 System.out.println(Y);
	 System.out.println();


	//Y=Q.transpose().times(Y);
	Y = Q.t() * Y;

	//	Y=Rs;
	//(Y - Rs).print("This is Y - RS");
	//Y.print("this is y dash: ");
	System.out.println("Y dash: ");
	 System.out.print(Y.transpose());
	 System.out.println();
//R.print();

	LSDTree* tr = new LSDTree(1000000, 32);
//	LSDTree tree = LSDTree(1000000, 32);
	//std::auto_ptr<LSDTree> tree(new LSDTree);

	tr->setRMatrix(R);

	tr->generateFirstlevel(Y.at(31, 0));

	//System.out.println(Y.get(31, 0));
	//tree.printTree();
	//tree.printcurrentlevel();
	//System.out.println("***************");

	for (int i = 0; i < 31; i++) {
		//	System.out.println("************************** After level : "+(i+1));
		//System.out.println("Y index: "+i+1);
		tr->generateNextlevel(Y.at(31 - (i + 1), 0), i + 1);

		//tree.printTree();
		//tree.printcurrentlevel();
		//	System.out.println("*************************************");
	}
	//System.out.println("*************************************");

	 minNode = tr->getMinnode();
	//cout << "***********************" <<minNode.getvalue()<< "  "<<minNode.getNode_S()<<endl;
	//System.out.println(minNode.getNode_S()+" : "+minNode.getvalue());
	for (int i = 0; i < 32; i++) {
		//	System.out.println(i+" "+minNode.getvalue());
		y[i] = minNode.getNode_S();

		minNode = (minNode.getparent());
	}

	delete tr;
	//	cout << "***********************" << endl;
*/

	CreateTree *ct;

	ct= new CreateTree();
//	cout << "***********************" << endl;
	y=ct->runtree(Y,H);
	delete ct;
	return y;
}

void Equalizer::RS(mat S) {
	/*QRDecomposition QRD = new QRDecomposition(H);
	 Q=QRD.getQ();
	 R=QRD.getR();*/

	arma::qr(Q, R, H);
	//Q.print("This is Q : ");
	//R.print(""This is R : ");
	Rs = R * S;
	cout << "This is RS: " << endl << Rs.t() << endl;

}
mat Equalizer::getRS(mat S) {
	/*QRDecomposition QRD = new QRDecomposition(H);
	 Q=QRD.getQ();
	 R=QRD.getR();*/

	arma::qr(Q, R, H);
	//Q.print("This is Q : ");
	//R.print(""This is R : ");
	Rs = R * S;
	return Rs;

}
mat Equalizer::getQtRS(mat S) {
	/*QRDecomposition QRD = new QRDecomposition(H);
	 Q=QRD.getQ();
	 R=QRD.getR();*/

	arma::qr(Q, R, H);
	//Q.print("This is Q : ");
	//R.print(""This is R : ");
	Rs = R * S;
	return (Q.t().i()) * Rs;

}
