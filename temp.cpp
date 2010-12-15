#include "mex.h"
#include "peopt.h"
#include "simple_matching.h"

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. oldY
    // 2. oldS
    // 3. p
    if(nrhs!=2)
	mexErrMsgTxt("The function BFGS requires 3 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function BFGS requires 1 output.");

    // Get the size of info 
    const int minfo=mxGetM(prhs[0]);
    const int ninfo=mxGetN(prhs[0]);

    // Get the size of our direction of evaluation, p
    const int mp=mxGetM(prhs[1]);
    const int np=mxGetN(prhs[1]);

    // Verify the sizes
    if(minfo!=mp)
    	mexErrMsgTxt("The matrix containing the gradient information and "
	 "old step information must contain the same number of rows as the "
	 "direction of evaluation.");
    if(ninfo%2 != 0)
     	mexErrMsgTxt("The matrix containing the gradient and old step "
	    "information must contain and even number of columns.");
    if(np!=1)
    	mexErrMsgTxt("The direction of evaluation must be a column vector.");

    // Get the data for this function.  Convert it to a vector.
    const double* info0=mxGetPr(prhs[0]);
    list < vector <double> > info;
    for(int j=0;j<ninfo;j++){
    	vector <double> tmp(info0+minfo*j,info0+minfo*(j+1));
    	info.push_back(tmp);
    }

    // Extract the direction we evaluate the Hessian 
    const double* p0=mxGetPr(prhs[1]);
    vector <double> p(p0,p0+np);
    
    // Allocate memory for the workspaces
    list < vector <double> > workU(0);

    // Allocate memory for the result of evaluation
    vector <double> result(mp);

    // Evaluate the BFGS preconditioner based on this data
    (Preconditioners::BFGS < vector <double> > (info))(p,result);

    // Create the output
    plhs[0] = mxCreateDoubleMatrix(mp,np,mxREAL);
    double* result0=mxGetPr(plhs[0]);
    for(int i=0;i<mp;i++)
    	result0[i]=result[i];
}
