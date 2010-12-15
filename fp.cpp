#include "mex.h"
#include "peopt.h"
#include "simple_matching.h"

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. i
    // 2. d
    // 3. y
    // 4. eta
    if(nrhs!=4)
	mexErrMsgTxt("The function fp requires 4 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function fp requires 1 output.");

    // Get the size of d
    const int m=mxGetM(prhs[1]);
    const int n=mxGetN(prhs[1]);

    // Get the size of y
    const int my=mxGetM(prhs[2]);
    const int ny=mxGetN(prhs[2]);
    
    // Get the size of eta 
    const int meta=mxGetM(prhs[2]);
    const int neta=mxGetN(prhs[2]);

    // Check that y has only 1 column
    if(ny!=1)
    	mexErrMsgTxt("The point of evaluation must have a single column.");
    // Check that eta has only 1 column
    if(ny!=1)
    	mexErrMsgTxt("The direction of differentiation must have a single "
	    "column.");
    if(my!=m) 
    	mexErrMsgTxt("The point of evaluation and the data must have the "
	    "same number of rows.");

    // Get which objective function we're working with
    const int i=(int)mxGetScalar(prhs[0]);

    // Make sure this index is valid.  It must correspond to a col of the data
    if(i<1 || i>n)
    	mexErrMsgTxt("The index i is invalid.");

    // Get the data for this function.  Convert it to a vector.
    const double* d=mxGetPr(prhs[1]);
    vector < vector <double> > dd;
    for(int j=0;j<n;j++){
    	vector <double> tmp(d+m*j,d+m*(j+1));
    	dd.push_back(tmp);
    }

    // Get the point we evaluate this function at.  Convert it to a vector.
    const double* y=mxGetPr(prhs[2]);
    vector <double> yy(y,y+mxGetM(prhs[2]));
    
    // Get the direction of differentiation.  Convert it to a vector.
    const double* eta=mxGetPr(prhs[3]);
    vector <double> eta_eta(eta,eta+mxGetM(prhs[3]));
	
    // Create the matching operator
    BasicMatching op(dd);

    // Evaluate the objective function at this point
    vector <double> result(m);
    op.p(i-1,yy,eta_eta,result);

    // Create the output from the function
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
    double *output=(double*)mxGetData(plhs[0]);
    for(int j=0;j<m;j++)
    	output[j]=result[j];
}
