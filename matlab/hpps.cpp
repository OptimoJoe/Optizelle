#include "mex.h"
#include "core.h"
#include "simple_matching.h"
#include<iostream>

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. A (cell array)
    // 2. B
    // 3. b (cell array)
    // 4. i
    // 5. u
    // 6. eta
    // 7. xi
    if(nrhs!=7)
	mexErrMsgTxt("The function hpps requires 7 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function hpps requires 1 output.");

    // Get the size of the collection of As
    const int mA=mxGetM(prhs[0]);
    const int nA=mxGetN(prhs[0]);

    // Get the size of B
    const int mB=mxGetM(prhs[1]);
    const int nB=mxGetN(prhs[1]);

    // Get the size of the collection of bs
    const int mb=mxGetM(prhs[2]);
    const int nb=mxGetN(prhs[2]);

    // Get the size of our index
    const int mi=mxGetM(prhs[3]);
    const int ni=mxGetN(prhs[3]);

    // Get the size of our control
    const int mu=mxGetM(prhs[4]);
    const int nu=mxGetN(prhs[4]);
    
    // Get the size of the direction of differentiation
    const int meta=mxGetM(prhs[5]);
    const int neta=mxGetN(prhs[5]);
    
    // Get the size of the second direction of differentiation
    const int mxi=mxGetM(prhs[6]);
    const int nxi=mxGetN(prhs[6]);

    // Verify the sizes
    if(mA!=1)
    	mexErrMsgTxt("The cell array containing the linear operators must "
	    "contain only a single row.");
    if(mb!=1)
    	mexErrMsgTxt("The cell array containing the data must contain only "
	    "a single row.");
    if(mu!=1)
	mexErrMsgTxt("The control must contain a single row.");
    if(nu!=nA)
    	mexErrMsgTxt("The size of the control must match the number of linear "
	    "operators.");
    if(mi!=1 || ni!=1)
    	mexErrMsgTxt("The index must be a scalar.");
    if(mB!=nB)
    	mexErrMsgTxt("The linear operators must be square.");
    if(neta!=nu || meta!=mu)
    	mexErrMsgTxt("The size of the direction of differenetiation must match "
	    "the control.");
    if(nxi!=1 || mxi!=nB)
    	mexErrMsgTxt("The size of the second direction of differenetiation "
	    "must be a column vector with size that matches the linear "
	    "operators.");

    // Extract the linear operators
    if(!mxIsCell(prhs[0]))
    	mexErrMsgTxt("The first argument must be a cell array containing the "
	    "linear operators for the problem.");
    vector < vector<double> > A;
    for(int i=0;i<nA;i++){
    	// Extract the cell
	mxArray* AAA=mxGetCell(prhs[0],i);

	// Find its size
	const int mAAA=mxGetM(AAA);
	const int nAAA=mxGetM(AAA);

	// Verify this size aligns with B
	if(mAAA!=mB || nAAA!=nB)
	    mexErrMsgTxt("The size of the linear operators must match");

	// Extract a double pointer from this matrix
    	const double* AA=mxGetPr(AAA);
	A.push_back ( vector <double> (AA,AA+mAAA*nAAA));
    }

    // Extract the data
    if(!mxIsCell(prhs[2]))
    	mexErrMsgTxt("The third argument must be a cell array containing the "
	    "data for the problem.");
    vector < vector<double> > b;
    for(int i=0;i<nb;i++){
    	// Extract the cell
	mxArray* bbb=mxGetCell(prhs[2],i);

	// Find its size
	const int mbbb=mxGetM(bbb);
	const int nbbb=mxGetN(bbb);

	// Verify this size aligns with B and has one column
	if(mbbb!=mB || nbbb!=1)
	    mexErrMsgTxt("The size of the data must align with B and have "
		"a single column.");

	// Extract a double pointer from this matrix
    	const double* bb=mxGetPr(bbb);
	vector <double> tmp(bb,bb+mbbb);
	b.push_back (tmp);
    }

    // Extract the linear operator that doesn't depend on the control
    const double* BB=mxGetPr(prhs[1]);
    vector <double> B(BB,BB+mB*nB);

    // Extract the index
    const int i=(int)mxGetScalar(prhs[3]);

    // Check that the index does not exceed the size of b
    if(i>nb)
	mexErrMsgTxt("The index exceeds the number of pieces of data");
    
    // Extract the control
    const double* uu=mxGetPr(prhs[4]);
    vector <double> u(uu,uu+nu);

    // Extract the direction of differentiation
    const double* eta_eta=mxGetPr(prhs[5]);
    vector <double> eta(eta_eta,eta_eta+neta);

    // Extract the second direction of differentiation
    const double* xi_xi=mxGetPr(prhs[6]);
    vector <double> xi(xi_xi,xi_xi+mxi);

    // Form the linear operator
    BasicOp op(A,B,b);

    // Evaluate the operator
    vector <double> y;
    op.pps(i-1,u,eta,xi,y);

    // Create the output
    plhs[0] = mxCreateDoubleMatrix(1,nu,mxREAL);
    double* yy=mxGetPr(plhs[0]);
    for(int i=0;i<nu;i++)
    	yy[i]=y[i];
}
