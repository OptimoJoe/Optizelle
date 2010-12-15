#include "mex.h"
#include "peopt.h"
#include "simple_matching.h"

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. i
    // 2. d
    // 3. y
    if(nrhs!=3)
	mexErrMsgTxt("The function f requires 3 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function f requires 1 output.");

    // Get the size of d
    const int m=mxGetM(prhs[1]);
    const int n=mxGetN(prhs[1]);

    // Get the size of y
    const int my=mxGetM(prhs[2]);
    const int ny=mxGetN(prhs[2]);

    // Check that y has only 1 column and that d and y have the same
    // number of rows
    if(ny!=1)
    	mexErrMsgTxt("The point of evaluation must have a single columns.");
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
	
    // Create the matching operator
    BasicMatching op(dd);

    // Evaluate the objective function at this point
    vector <double> result(m);
    op(i-1,yy,result);

    // Create the output from the function
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
    double *output=(double*)mxGetData(plhs[0]);
    for(int j=0;j<m;j++)
    	output[j]=result[j];
}
