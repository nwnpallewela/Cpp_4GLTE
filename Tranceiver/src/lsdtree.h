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
#include <list>
#include <vector>

using namespace std;
using namespace arma;

bool rule(const std::vector<double> a, const std::vector<double> b) {
	return (a.at(0)) > (b.at(0));
}

class Node {
private:
	int S_values[8];
	int S;
	double value;
	Node *parent;
	vector<Node*> nextlevel;

// constructor
public:

	void setvalue(double value);
	void copyNode1(int S, double value);
	void copyNode1(int S, double value, Node *Parent);
	void copyNode(int S, double value, Node *Parent, vector<Node*> next);
	Node* findNode_with_S(int s);
	Node getNode_with_S(int s);
	double getvalue();
	int getNode_S();
	Node getparent();
	vector<Node*> getnextlevel();
	void addNodetonextlevel(int i, Node *parent, double value);
	void printnum();
	void printNumValue();
	void printvalue();
	void printParents();
	void printNodes();
	Node* getparent_pointer();
	Node* getchildS_pointer(int s);
	Node() {
		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}

		this->S = 0;
		this->parent = NULL;
		this->nextlevel = vector<Node*>();
		this->value = 0;
	}
	Node(int i) {
		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}

		this->S = S_values[i];
		this->parent = NULL;
		this->nextlevel = vector<Node*>();
		this->value = 0;
	}
	Node(int i, Node *parent) {
		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}

		this->S = S_values[i];
		this->parent = parent;
		this->nextlevel = vector<Node*>();
		this->value = 0;
	}

	Node(int i, Node *parent, double value) {
		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}
		this->S = S_values[i];
		this->parent = parent;
		this->nextlevel = vector<Node*>();
		this->value = value;
	}
	Node(int i, double value) {
		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}
		this->S = S_values[i];
		this->parent = NULL;
		this->nextlevel = vector<Node*>();
		this->value = value;
	}

// //////////////////////////////////////////////////

};

void Node::copyNode1(int S, double value, Node *Parent) {
	this->S = S;
	this->value = value;
	this->parent = Parent;
}
void Node::copyNode1(int S, double value) {
	this->S = S;
	this->value = value;
	this->parent = NULL;
}

void Node::copyNode(int S, double value, Node *Parent, vector<Node*> next) {
	//cout << "Point of problem copy node" << endl;
	this->S = S;
	this->value = value;
	this->parent = Parent;
	this->nextlevel = next;
}

void Node::setvalue(double value) {
	this->value = value;
}

Node* Node::findNode_with_S(int s) {
	Node *temp;
	for (int i = 0; i < nextlevel.size(); ++i) {

		if (nextlevel[i]->getNode_S() == s) {

			return nextlevel[i];
			//return temp;
		}

	}

	return temp;

}
Node Node::getNode_with_S(int s) {
	Node temp;
	for (int i = 0; i < nextlevel.size(); ++i) {

		if (nextlevel[i]->getNode_S() == s) {

			return *nextlevel[i];
			//return temp;
		}

	}

	return temp;

}

// //////////////////////////////////////////////////
double Node::getvalue() {
	return value;
}

int Node::getNode_S() {
	//cout << "Point of problem" << endl;
	return S;
}

Node Node::getparent() {
	return *parent;
}

Node* Node::getparent_pointer() {
	return parent;
}

vector<Node*> Node::getnextlevel() {
	return nextlevel;
}

void Node::addNodetonextlevel(int i, Node *parent, double value) {
	this->nextlevel.push_back((new Node(i, parent, value))); //( Node(i, parent, value));
}

// //////////////////////////////////////////////////
void Node::printnum() {
	cout << S << " : " << endl;

}

void Node::printvalue() {
	cout << "value : " << value << endl;

}

void Node::printNumValue() {
	cout << S << " : " << value << endl;
}

void Node::printParents() {
	Node temp = Node(1);
	temp = *parent;
	while (temp.value != 0) {
		cout << "(" << temp.getNode_S() << " : " << temp.getvalue() << ") "
				<< endl;
		temp = temp.getparent();
	}
}

void Node::printNodes() {

	if (nextlevel.empty()) {
		cout << "(" << S << " : " << value << ") ! " << endl;
	} else {

		for (std::vector<Node*>::iterator i = nextlevel.begin();
				i != nextlevel.end(); ++i) {
			cout << "  ---   ";
			//i.printNodes();
			(*(i.operator *())).printNodes();

		}
	}

	cout << "-->(" << S << " : " << value << ") " << endl;
}

Node* Node::getchildS_pointer(int s) {
	Node *temp;
	for (std::vector<Node*>::iterator i = nextlevel.begin();
			i != nextlevel.end(); ++i) {
		if ((*(i.operator *())).getNode_S() == s) {
			temp = (i.operator *());
		}

	}
	return temp;
}

////////////////////////////////////////////////////////////////////////////////////////////

class LSDTree {

private:

	Node root;
	Node minNode;
	int S_values[8];

	double error[64];
	double R; // Sphere r value
	int K;
	int precision; // upto which number of decimal places

	vector<vector<double> > level_values_next;
	vector<Node*> currentlevel;
	vector<Node*> nextlevel;
	mat R_mat;

	bool testSphere(int i, double r, double im);
	double calcValue(int i, double y, int level, Node node);

public:
	void printTree();
	void generateFirstlevel(double y);
	void setRMatrix(mat R);
	void generateNextlevel(double y, int level);
	Node getMinnode();
	void printcurrentlevel();

	LSDTree(double R, int K) {

		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}
		this->root = Node(0);
		this->root.setvalue(0);
		this->minNode = Node(0, 9999999999999.9);
		currentlevel = vector<Node*>();
		nextlevel = vector<Node*>();

		//level_values_next =  list<>();
		//this->error =  double[64];
		this->K = K;
		this->R = R;
		this->precision = 1000000;
	}
	LSDTree() {

		for (int var = 0; var < 8; ++var) {
			S_values[var] = 2 * var - 7;
		}
		this->root = Node(0);
		this->root.setvalue(0);
		this->minNode = Node(0, 9999999999999.9);
		currentlevel = vector<Node*>();
		nextlevel = vector<Node*>();

		//level_values_next =  list<>();
		//this->error =  double[64];
		this->K = 32;
		this->R = 1000000.0;
		this->precision = 1000000;
	}

};

void LSDTree::printTree() {

	root.printNodes();
}

void LSDTree::generateFirstlevel(double y) {

	//cout << "this is y " << y << endl;

	int count = 0;
	double temp2 = 0;
	int index = 0;
	vector<double> index_val;
	level_values_next.clear();
	/*
	 * currentlevel.clear(); nextlevel.clear();
	 */
	for (int i = 0; i < 8; i++) {
		vector<double> temp1(3);

		temp1[0] = (calcValue(i, y, 0, root));
		//	cout << i << " -- " << temp1[0] << endl;
		temp1[1] = (double) i;
		temp1[2] = (double) 0;
		level_values_next.push_back(temp1);
		count++;

	}

// Collections.sort(level_values);

	/*	Collections.sort(level_values_next, new Comparator<Double[]>() {

	 @Override
	 public int compare(Double[] arg0, Double[] arg1) {
	 // TODO Auto-generated method stubnt

	 if (Double.compare(arg0[0],arg1[0])>0) {
	 return -1;
	 } else {
	 return 1;
	 }

	 }
	 });
	 */

	std::sort(level_values_next.begin(), level_values_next.end(), rule);

	//cout<< S_values[(int)(level_values_next[level_values_next.size() - 1].at(1))] <<" - "<<level_values_next[level_values_next.size() - 1].at(0)<<endl;

	for (int i = level_values_next.size() - 1; i >= 0; i--) {
		if (i >= (count - K)) {

			index_val = level_values_next.at(i);
			index = (int) index_val[1];
			temp2 = index_val[0];
			// System.out.println(temp);
			// System.out.println(index);
			// System.out.println(temp);

			// System.out.println(node_index);

			root.addNodetonextlevel(index, &root, temp2);

		}

	}

	int count_cl = 0;

	/*	for (Iterator<Node> iterator = root.getnextlevel().iterator();
	 iterator.hasNext();) {
	 Node node = (Node) iterator.next();
	 currentlevel.add(count_cl, new Node(0, null));
	 currentlevel.get(count_cl).copyNode(node.getNode_S(), node.getvalue(),
	 node.getparent(), node.getnextlevel());
	 count_cl++;
	 }*/
//	cout<<root.getnextlevel().size()<<endl;
	/*	for (std::vector<Node*>::iterator i = root.getnextlevel().begin();
	 i != root.getnextlevel().end(); ++i) {

	 cout<< (i.operator *())->getNode_S()<<endl;
	 Node *temp = new Node(0);
	 cout<<"Point of problem for loop"<<endl;
	 cout<<((*i.operator *()).getNode_S())<<endl;
	 cout<<"Point of problem for loop"<<endl;
	 cout<<(*i.operator *()).getvalue()<<endl;
	 cout<<(*(*i.operator *()).getparent_pointer()).getNode_S()<<endl;
	 cout<<(*i.operator *()).getnextlevel().size()<<endl;


	 (*temp).copyNode((*i.operator *()).getNode_S(), (*i.operator *()).getvalue(),
	 ((*i.operator *()).getparent_pointer()),
	 (*i.operator *()).getnextlevel());



	 currentlevel.push_back(temp);
	 count_cl++;

	 }*/
	// std::vector<int>::size_type sz = myvector.size();
	for (unsigned var = 0; var < root.getnextlevel().size(); ++var) {
		Node *temp = new Node(0);
		//	root.getnextlevel()[var]->printNumValue();
		(*temp).copyNode(root.getnextlevel()[var]->getNode_S(),
				root.getnextlevel()[var]->getvalue(),
				(root.getnextlevel()[var]->getparent_pointer()),
				root.getnextlevel()[var]->getnextlevel());
		//temp->printNumValue();
		currentlevel.push_back(temp);
		//	cout << currentlevel.size() << endl;
		count_cl++;
	}

	/*
	 * System.out.println("************************");
	 * System.out.println("Y : " + y); for (Iterator<Node> iterator =
	 * currentlevel.iterator(); iterator .hasNext();) { Node node = (Node)
	 * iterator.next(); node.printNumValue(); }
	 */

	/*
	 * for (Iterator<Double[]> iterator = level_values_next.iterator();
	 * iterator .hasNext();) { Double[] val = iterator.next();
	 * System.out.println(val[0] + " : " + val[1] + " " + val[2]); }
	 */

}

void LSDTree::setRMatrix(mat R) {
	R_mat = R;
}

void LSDTree::generateNextlevel(double y, int level) {

	int count = 0;
	double temp2 = 0;
	int index = 0;
	int node_index = 0;
	vector<double> index_val;
	double current_node_val = 0;

	level_values_next.clear();

	/*
	 * for (int i = 0; i < 8; i++) { error[i] = calcValue(i, y, level);
	 *
	 * System.out.println("error = " +real[i]+" + "+complex[i]+"i -- "+in
	 * [0]+" + "+in[1]+"i : "+error[i]);
	 *
	 * }
	 */

	/*
	 * System.out.println("************************" + level);
	 * System.out.println("Y : " + y);
	 */

	for (int j = 0; j < currentlevel.size(); j++) {
		//	cout << "*******************" << currentlevel.size() << " " << endl;

		current_node_val = currentlevel.at(j)->getvalue();

		for (int i = 0; i < 8; i++) {

			//	System.out.print(" sent node : ");
			/*
			 * currentlevel.get(j).getparent()
			 * .findNode_with_S(currentlevel.get(j).getNode_S())
			 * vector<.printNumValue();
			 */
			vector<double> temp(3);
			/*		(calcValue(i, y, level,
			 currentlevel.at(j).getparent().findNode_with_S(
			 currentlevel.at(j).getNode_S()))
			 + current_node_val), (double) i, (double) j };

			 */
			//vector<double> temp1 ;
			//cout << "    " << currentlevel[j]->getNode_S() << "   ";
			temp[0] = (calcValue(i, y, level,
					currentlevel[j]->getparent().getNode_with_S(
							currentlevel[j]->getNode_S())) + current_node_val);
			temp[1] = (double) i;
			temp[2] = (double) j;
			level_values_next.push_back(temp);
			count++;
			//cout << "*******************" << endl;
		}

	}

//	System.out.println("level_values_next"+level_values_next.size());

	/*	Collections.sort(level_values_next, new Comparator<Double[]>() {

	 @Override
	 public int compare(Double[] arg0, Double[] arg1) {
	 // TODO Auto-generated method stubnt

	 return arg1[0].compareTo(arg0[0]);

	 }
	 });*/

	std::sort(level_values_next.begin(), level_values_next.end(), rule);
//	cout<< S_values[(int)(level_values_next[level_values_next.size() - 1].at(1))] <<" - "<<level_values_next[level_values_next.size() - 1].at(0)<<endl;
	//cout<<"best "<<level_values_next[0].at(0)<<endl;

	/*
	 * for (Iterator<Double[]> iterator = level_values_next.iterator();
	 * iterator .hasNext();) { Double[] val = iterator.next();
	 * System.out.println(val[0]+" : "+val[1]+" "+val[2]); }
	 */
	int count_k = 0;
//	cout << "this is values" << endl;

	for (int i = level_values_next.size() - 1; i >= 0; --i) {
		if (count_k < (K)) {
			//	cout << level_values_next[i].at(0) << endl;
			count_k++;
			index_val = level_values_next.at(i);

			index = (int) index_val[1];
			temp2 = index_val[0];
			node_index = (int) index_val[2];
			// System.out.println(index+" "+node_index+" "+temp);
			// System.out.println(node_index);
			/*
			 * currentlevel.get(node_index).addNodetonextlevel(index,
			 * currentlevel.get(node_index), temp);
			 */
			/*
			 cout<<(*currentlevel[node_index]->getparent().findNode_with_S(
			 currentlevel[node_index]->getNode_S())).getNode_S()<<" "<<(*currentlevel[node_index]->getparent().findNode_with_S(
			 currentlevel[node_index]->getNode_S())).getparent().getNode_S()<<endl<<endl;
			 */

			(*currentlevel[node_index]->getparent().findNode_with_S(
					currentlevel[node_index]->getNode_S())).addNodetonextlevel(
					index,
					currentlevel[node_index]->getparent().findNode_with_S(
							currentlevel[node_index]->getNode_S()), temp2);

			/*cout<<"next level size: "<<(*currentlevel[node_index]->getparent().findNode_with_S(
			 currentlevel[node_index]->getNode_S())).getNode_S()<<" "<<(*currentlevel[node_index]->getparent().findNode_with_S(
			 currentlevel[node_index]->getNode_S())).getnextlevel().size()<<endl<<endl;
			 */
			/*
			 (*(*(currentlevel[node_index]->getparent_pointer())).findNode_with_S(
			 currentlevel[node_index]->getNode_S())).addNodetonextlevel(
			 index,
			 ((*(currentlevel[node_index]->getparent_pointer())).findNode_with_S(
			 currentlevel[node_index]->getNode_S())),
			 temp2);*/
			// System.out.println(count_k+" : "+temp);
		}

	}
	nextlevel.clear();
	minNode.copyNode1(0, 99999999.9);			//null
// minNode.copyNode(node1.getNode_S(), node1.getvalue(),
// node1.getparent());
// minNode.setvalue(999999999999999.9);

	int count_nl = 0;
	/*for (Iterator<Node> iterator = currentlevel.iterator(); iterator.hasNext();
	 ) {
	 Node node = (Node) iterator.next();
	 for (Iterator<Node> iterator2 = node.getnextlevel().iterator();
	 iterator2.hasNext();) {
	 Node node1 = (Node) iterator2.next();
	 if (minNode.getvalue() > node1.getvalue()) {
	 minNode.copyNode(node1.getNode_S(), node1.getvalue(),
	 node1.getparent());
	 }
	 // nextlevel.add(node1);
	 nextlevel.push_back(  Node(0, NULL));
	 nextlevel.pop_back().copyNode(node1.getNode_S(),
	 node1.getvalue(), node1.getparent(), node1.getnextlevel());
	 count_nl++;

	 }

	 }*/

	/*	for (std::vector<Node*>::iterator i = currentlevel.begin();
	 i != currentlevel.end(); ++i) {

	 for (std::vector<Node*>::iterator i1 =
	 (*i.operator *()).getnextlevel().begin();
	 i1 != (*i.operator *()).getnextlevel().end(); ++i1) {
	 if (minNode.getvalue() > (*i1.operator *()).getvalue()) {

	 minNode.copyNode1((*i1.operator *()).getNode_S(),
	 (*i1.operator *()).getvalue(),
	 ((*i1.operator *()).getparent_pointer()));
	 }
	 // nextlevel.add(node1);

	 Node *temp = new Node(0);

	 (*temp).copyNode((*i1.operator *()).getNode_S(),
	 (*i1.operator *()).getvalue(),
	 ((*i1.operator *()).getparent_pointer()),
	 (*i1.operator *()).getnextlevel());

	 nextlevel.push_back(temp);
	 count_nl++;

	 }
	 }*/
	//cout << "current level size" << currentlevel.size() << endl;
	for (unsigned i = 0; i < currentlevel.size(); ++i) {
		/*cout << "current level next size " << i << " "
		 << currentlevel[i]->getparent().getNode_with_S(currentlevel[i]->getNode_S()).getnextlevel().size() << endl;
		 */for (unsigned j = 0;
				j
						< currentlevel[i]->getparent().getNode_with_S(
								currentlevel[i]->getNode_S()).getnextlevel().size();
				++j) {

			if (minNode.getvalue()
					> currentlevel[i]->getparent().getNode_with_S(
							currentlevel[i]->getNode_S()).getnextlevel()[j]->getvalue()) {

				minNode.copyNode1(
						currentlevel[i]->getparent().getNode_with_S(
								currentlevel[i]->getNode_S()).getnextlevel()[j]->getNode_S(),
						currentlevel[i]->getparent().getNode_with_S(
								currentlevel[i]->getNode_S()).getnextlevel()[j]->getvalue(),
						currentlevel[i]->getparent().getNode_with_S(
								currentlevel[i]->getNode_S()).getnextlevel()[j]->getparent_pointer());
			}
			// nextlevel.add(node1);

			Node *temp = new Node(0);

			(*temp).copyNode(
					currentlevel[i]->getparent().getNode_with_S(
							currentlevel[i]->getNode_S()).getnextlevel()[j]->getNode_S(),
					currentlevel[i]->getparent().getNode_with_S(
							currentlevel[i]->getNode_S()).getnextlevel()[j]->getvalue(),
					currentlevel[i]->getparent().getNode_with_S(
							currentlevel[i]->getNode_S()).getnextlevel()[j]->getparent_pointer(),
					currentlevel[i]->getparent().getNode_with_S(
							currentlevel[i]->getNode_S()).getnextlevel()[j]->getnextlevel());

			nextlevel.push_back(temp);
			count_nl++;

		}

	}
	//cout << "next level size" << nextlevel.size() << " " << count_nl << endl;
// Node temp_val=new Node(0, null);
// temp_val=minNode;
// System.out.print(" this level best  value : ");
// System.out.println(minNode.getvalue());
	currentlevel.clear();
	int count_cl = 0;

	/*for (Iterator<Node> iterator = nextlevel.iterator(); iterator.hasNext();) {
	 Node node = (Node) iterator.next();
	 currentlevel.add(count_cl, new Node(0, null));
	 currentlevel.get(count_cl).copyNode(node.getNode_S(), node.getvalue(),
	 node.getparent(), node.getnextlevel());

	 count_cl++;
	 }*/
	/*for (std::vector<Node*>::iterator i = nextlevel.begin();
	 i != nextlevel.end(); ++i) {
	 //currentlevel.push_back(Node(0));

	 Node *temp = new Node(0);

	 (*temp).copyNode((*i.operator *()).getNode_S(),
	 (*i.operator *()).getvalue(),
	 ((*i.operator *()).getparent_pointer()),
	 (*i.operator *()).getnextlevel());
	 currentlevel.push_back(temp);
	 count_cl++;
	 }*/
	//
	//cout<<"next level size "<<nextlevel.size()<<endl;
	for (unsigned var = 0; var < nextlevel.size(); ++var) {
		Node *temp = new Node(0);
		//	root.getnextlevel()[var]->printNumValue();
		(*temp).copyNode(nextlevel[var]->getNode_S(),
				nextlevel[var]->getvalue(),
				(nextlevel[var]->getparent_pointer()),
				nextlevel[var]->getnextlevel());

		currentlevel.push_back(temp);
		count_cl++;
	}

// int count_level = 1;
// System.out.println("**************");
	/*
	 * for (Iterator<Node> iterator = currentlevel.iterator(); iterator
	 * .hasNext();) { Node node = (Node) iterator.next();
	 * System.out.print(count_level+" : ");
	 * //node.getparent().printNumValue(); node.printNumValue();
	 * count_level++; }
	 */

}

Node LSDTree::getMinnode() {
	return minNode;
}
void LSDTree::printcurrentlevel() {
	int count = 1;
	/*for (Iterator<Node> iterator = currentlevel.iterator(); iterator.hasNext();
	 ) {
	 Node node = (Node) iterator.next();

	 }*/

	for (std::vector<Node*>::iterator i = currentlevel.begin();
			i != currentlevel.end(); ++i) {
		cout << count << " : ";
		(*i.operator *()).printNumValue();
		count++;

	}

}

double LSDTree::calcValue(int i, double y, int level, Node node) {

	double value = 0;
	double sum = 0.0;
	int S_temp = S_values[i];
	Node temp = Node(1);
	temp.copyNode1(node.getNode_S(), node.getvalue(), node.getparent_pointer());
//	cout << temp.getNode_S() << endl;
	if (level == 0) {

		value = (y - S_values[i] * R_mat(31, 31));
		/*cout << "calc value : " << y << " - " << S_values[i] << " * "
		 << R_mat(31, 31) << "(" << S_values[i] * R_mat(31, 31) << ") = "
		 << value << endl;*/
		value = value * value;

	} else {
		//System.out.print("  level " + (level) + " : ");

		sum = S_temp * R_mat(31 - level, 31 - level);

		//System.out.println("S_value " + S_temp + "\t j:-1");
		//	System.out.print("          : ");
		S_temp = temp.getNode_S();
		for (int j = 0; j < level; j++) {
			sum = S_temp * R_mat.at(31 - level, 31 - level + (j + 1)) + sum;

			//	System.out.println("S_value " + S_temp + "\t j: " + j + " \t "
			//		+ temp.getNode_S() + " " + temp.getvalue());

			//if (j != 0 && j < level) {

			if (j < level) {

				//	cout<<temp.getparent().getNode_S()<<endl;
				//	cout << "*******************" << endl;
				temp.copyNode1(temp.getparent().getNode_S(),
						temp.getparent().getvalue(),
						temp.getparent().getparent_pointer());

				S_temp = temp.getNode_S();
			}

			/* else if (j == 0) {

			 S_temp = temp.getNode_S();
			 }*/
			//	System.out.print("          : ");
		}
		//	System.out.println();
		value = (y - sum) * (y - sum);
	}

	return (value);
}

bool LSDTree::testSphere(int i, double r, double im) {

	return true;

}
class CreateTree {
private:
	LSDTree tree;
	mat H;
	mat Q;
	mat R;
	mat Rs;
	Node minNode;
	double y[32];

public:
	CreateTree() {


		tree = LSDTree(1000000, 32);
		H = mat(32, 32, fill::zeros);
		Q = mat(32, 32, fill::zeros);
		R = mat(32, 32, fill::zeros);
		Rs = mat(32, 1, fill::zeros);

		minNode = Node();

	//	double y[32];
	}
	~CreateTree() {



		//	double y[32];
		}
	double* runtree(mat Y,mat Hin);

};

double* CreateTree::runtree(mat Y,mat Hin) {
H=Hin;
	arma::qr(Q, R, H);

	Y = Q.t() * Y;
	//cout<<"***********"<<endl;
	tree.setRMatrix(R);

	tree.generateFirstlevel(Y.at(31, 0));

	//System.out.println(Y.get(31, 0));
	//tree.printTree();
	//tree.printcurrentlevel();
	//System.out.println("***************");

	for (int i = 0; i < 31; i++) {
		//	System.out.println("************************** After level : "+(i+1));
		//System.out.println("Y index: "+i+1);
		tree.generateNextlevel(Y.at(31 - (i + 1), 0), i + 1);

		//tree.printTree();
		//tree.printcurrentlevel();
		//	System.out.println("*************************************");
	}
	//System.out.println("*************************************");

	minNode = tree.getMinnode();
	//cout << "***********************" <<minNode.getvalue()<< "  "<<minNode.getNode_S()<<endl;
	//System.out.println(minNode.getNode_S()+" : "+minNode.getvalue());
	for (int i = 0; i < 32; i++) {
		//	System.out.println(i+" "+minNode.getvalue());
		y[i] = minNode.getNode_S();

		minNode = (minNode.getparent());
	}

	//delete tree;
	//	cout << "***********************" << endl;

	return y;
}

