#ifndef SIMPLE_MATCHING_H
#define SIMPLE_MATCHING_H

#include <vector>
extern "C"{
    #include <clapack.h>
}

#include "diff_operator.h"

using namespace std;

// This implements the function f_i(y)=y-d_i
class BasicMatching : public DiffOperator < vector <double>, vector <double> > {
public:
    vector < vector <double> > d;	
    BasicMatching(vector < vector <double> >& d_) : d(d_) {};

    // Basic application
    void operator () (const int i,const vector <double>& x,vector <double>& y)
    	const;

    // Derivative of the operator
    void p(const int i,const vector <double>& x,const vector <double>& eta,
	vector<double>& y) const;
    
    // Derivative of the adjoint 
    void ps(const int i,const vector <double>& x,const vector <double>& xi,
	vector<double>& y) const;
    
    // Second derivative of the adjoint 
    void pps(const int i,
	const vector <double>& x,
	const vector <double>& eta,
	const vector <double>& xi,
	vector<double>& result) const;

    // Maximum index for the above functions
    int max_index() const{return d.size();}
};

// This implements the solution operator for u1*A1+u2*A2+...+umAm+B=b_i
class BasicOp : public DiffOperator < vector <double>, vector <double> > {
public:
    vector < vector <double> > A;
    vector <double> B;
    vector < vector <double> > b; 
    
    BasicOp(
	vector < vector <double> >& A_,
	vector < double >& B_,
	vector < vector <double> >& b_)
	: A(A_),B(B_),b(b_) {};
    
    // Computes the operator eta1A1+...+etamAm 
    void getAeta(const vector <double>& u,vector <double>& AA) const;

    // Computes the operator u1A1+...+umAm+B
    void getForwardOp(const vector <double>& u,vector <double>& AA) const;
    
    // Basic application
    void operator () (const int i,const vector <double>& u,vector <double>& y)
    	const;
    
    // Derivative of the operator 
    void p(const int i,const vector <double>& u,const vector <double>& eta,
	vector<double>& y) const;
    
    // Derivative of the adjoint 
    void ps(const int i,const vector <double>& u,const vector <double>& xi,
	vector<double>& y) const;
    
    // Second derivative of the adjoint 
    void pps(const int i,
	const vector <double>& u,
	const vector <double>& eta,
	const vector <double>& xi,
	vector<double>& result) const; 

    // Maximum index for the above functions
    int max_index() const {return b.size();}
};
#endif
