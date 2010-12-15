#include <vector>
extern "C"{
    #include <clapack.h>
}
#include "simple_matching.h"

using namespace std;

namespace Operations{
    // Scalar multiple
    template <>
    void scal (const double alpha,vector <double>& x){
    	cblas_dscal(x.size(),alpha,&(const_cast< vector<double>& >(x)[0]),1);
    }

    // Copy
    template <>
    void copy (const vector <double>& from, vector <double>& to){
    	to=from;
    }

    // alpha*x+y
    template <>
    void axpy(
	const double alpha,
	const vector <double>& x,
	vector <double>& y
    ){
    	cblas_daxpy(x.size(),alpha,&(x[0]),1,&(y[0]),1);
    }

    // Takes some element and sets it to all zeros
    template <>
    void zero(vector <double>& x){
    	scal(0.,x);
    }

    // Inner product
    template <>
    double innr(const vector <double>& x,const vector <double>& y){
    	return cblas_ddot(x.size(),&(x[0]),1,&(y[0]),1);
    }
}

// Basic application
void BasicMatching::operator () (
    const int i,
    const vector <double>& x,
    vector <double>& y
) const{
    y=x;
    Operations::axpy(-1.,d[i],y);
}

// Derivative of the operator
void BasicMatching::p(
    const int i,
    const vector <double>& x,
    const vector <double>& eta,
    vector<double>& y
) const{
    y=eta;
}
    
// Derivative of the adjoint 
void BasicMatching::ps(
    const int i,
    const vector <double>& x,
    const vector <double>& xi,
    vector<double>& y
) const{
    y=xi;
}
    
// Second derivative of the adjoint 
void BasicMatching::pps(const int i,
    const vector <double>& x,
    const vector <double>& eta,
    const vector <double>& xi,
    vector<double>& result 
) const{
    Operations::zero(result);
}

    
// Computes the operator eta1A1+...+etamAm 
void BasicOp::getAeta(const vector <double>& u,vector <double>& AA) const{
    // Determine the size of u
    int m=u.size();

    // Multply the Ai by Ui
    for(int j=0;j<m;j++)
	Operations::axpy(u[j],A[j],AA);
}

// Computes the operator u1A1+...+umAm+B
void BasicOp::getForwardOp(const vector <double>& u,vector <double>& AA) const{
    // Determine the size of u
    int m=u.size();

    // Multiply the Ai by Ui
    getAeta(u,AA);

    // Add in B
    Operations::axpy(1.,B,AA);
}
    
// Basic application
void BasicOp::operator () (
    const int i,
    const vector <double>& u,
    vector <double>& y
) const{
    // Determine the size of each A squared
    int nsq=A[0].size();

    // Get the whole operator 
    vector <double> AA(nsq,0.);
    getForwardOp(u,AA);

    // Determine the size of the rhs (one side of each A)
    int m=b[0].size();

    // Solve the linear system AA y = b_i
    y=b[i];
    int ipiv[m];
    int info=clapack_dgesv(CblasColMajor,m,1,&(AA[0]),m,
	ipiv,&(y[0]),m);
}
    
// Derivative of the operator 
void BasicOp::p(
    const int i,
    const vector <double>& u,
    const vector <double>& eta,
    vector<double>& y
) const{

    // Get the size of each Ai squared
    int nsq=A[0].size();

    // Compile the whole operator
    vector <double> AA(nsq,0.);
    getForwardOp(u,AA);

    // Compile A(eta)
    vector <double> Aeta(nsq,0.);
    getAeta(eta,Aeta);

    // Get the size of each right hand size (size of the state)
    int m=b[0].size();

    // Compute the solution operator
    vector <double> sol(m);
    (*this)(i,u,sol);

    // Apply -A(eta) to the solution operator.  Store the result in y
    y.resize(m);
    cblas_dgemv(CblasColMajor,CblasNoTrans,m,m,-1.,&(Aeta[0]),m,
	&(sol[0]),1,0.,&(y[0]),1);

    // Solve the linear system where the rhs is the compute y from above
    int ipiv[m];
    int info=clapack_dgesv(CblasColMajor,m,1,&(AA[0]),m,
	ipiv,&(y[0]),m);
}
    
// Derivative of the adjoint 
void BasicOp::ps(
    const int i,
    const vector <double>& u,
    const vector <double>& xi,
    vector<double>& y
) const{
    // Get the size of each Ai squared
    int nsq=A[0].size();

    // Compile the whole operator
    vector <double> AA(nsq,0.);
    getForwardOp(u,AA);

    // Determine the size of the state variable	
    int m=xi.size();

    // Determine the size of the control variable
    int n=u.size();

    // Solve an adjoint problem where the rhs is given by xi (the direction
    // of interest.)  Store this in nu.
    vector <double> nu=xi;
    int ipiv[m];
    int info=clapack_dgesv(CblasRowMajor,m,1,&(AA[0]),m,
	ipiv,&(nu[0]),m);
    
    // Apply the operator -(A(.)h(u))^* to nu 
    vector <double> sol(m);
    vector <double> tmp(m);
    y.resize(n);

    // Compute the solution operator
    (*this)(i,u,sol);

    // Compute the action of the derivative's adjoint
    for(int j=0;j<n;j++){

	// Apply Aj to the solution.  Store it in tmp.
	Operations::zero(tmp);
	cblas_dgemv(CblasColMajor,CblasNoTrans,m,m,1.,&(A[j][0]),m,
	    &(sol[0]),1,1.,&(tmp[0]),1);

	// Find the negative of the inner product between tmp and nu
	y[j]=-Operations::innr(tmp,nu);
    }
}
    
// Second derivative of the adjoint 
void BasicOp::pps(
    const int i,
    const vector <double>& u,
    const vector <double>& eta,
    const vector <double>& xi,
    vector<double>& result 
) const{
    // Get the size of each Ai squared
    int nsq=A[0].size();

    // Compile the whole operator
    vector <double> AA(nsq,0.);
    getForwardOp(u,AA);

    // Determine the size of the state variable	
    int m=xi.size();
    
    // Determine the size of the control variable
    int n=u.size();
    
    // Solve an adjoint problem where the rhs is given by xi (the direction
    // of interest.)  Store this in nu.
    vector <double> nu=xi;
    int ipiv[m];
    int info=clapack_dgesv(CblasRowMajor,m,1,&(AA[0]),m,
	ipiv,&(nu[0]),m);
    
    // Apply the operator -(A(.)h'(u)eta)^* to nu.  Store the result in
    // tmp3;
    vector <double> tmp(m);
    vector <double> tmp2(m);
    vector <double> tmp3(n);
    for(int j=0;j<n;j++){
	// Compute h'(u)eta.  Store the result in tmp
	p(i,u,eta,tmp);

	// Apply Aj to tmp.  Store it in tmp2.
	Operations::zero(tmp2);
	cblas_dgemv(CblasColMajor,CblasNoTrans,m,m,1.,&(A[j][0]),m,
	    &(tmp[0]),1,1.,&(tmp2[0]),1);

	// Find the inner product between tmp2 and nu
	tmp3[j]=-Operations::innr(tmp2,nu);
    }

    // Apply the operator -A(eta)* to nu.  Store the result in tmp;
    vector <double> Aeta(nsq);
    getAeta(eta,Aeta);
    cblas_dgemv(CblasColMajor,CblasTrans,m,m,-1.,&(Aeta[0]),m,
	&(nu[0]),1,0.,&(tmp[0]),1);

    // Apply the operator h'(u)^* to tmp.  Store the result in result
    ps(i,u,tmp,result);

    // Add together result and tmp3.  Store the result in result.
    Operations::axpy(1.,tmp3,result);	
}

void pe_error(char message[]){
    printf("%s\n",message);
}
