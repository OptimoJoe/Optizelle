#include "mex.h"
#include "peopt.h"
#include "simple_matching.h"

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. oldY
    // 2. oldS
    // 3. p
    if(nrhs!=3)
	mexErrMsgTxt("The function BFGSinv requires 3 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function BFGSinv requires 1 output.");

    // Get the size of the old delta gradient information 
    const int moldY=mxGetM(prhs[0]);
    const int noldY=mxGetN(prhs[0]);

    // Get the size of the old delta step information 
    const int moldS=mxGetM(prhs[1]);
    const int noldS=mxGetN(prhs[1]);

    // Get the size of our direction of evaluation, p
    const int mp=mxGetM(prhs[2]);
    const int np=mxGetN(prhs[2]);

    // Verify the sizes
    if(moldY!=mp)
    	mexErrMsgTxt("The matrix containing the gradient information and "
	 "old step information must contain the same number of rows as the "
	 "direction of evaluation.");
    if(moldY != moldS && noldY !=noldS)
     	mexErrMsgTxt("The matrix containing the old gradient deltas and "
	    "the matrix containing the old step deltas must be the same "
	    "size.");
    if(np!=1)
    	mexErrMsgTxt("The direction of evaluation must be a column vector.");

    // Get the old gradient deltas for this function.  Convert it to a vector.
    const double* y=mxGetPr(prhs[0]);
    list < vector <double> > oldY;
    for(int j=0;j<noldY;j++){
    	vector <double> tmp(y+moldY*j,y+moldY*(j+1));
    	oldY.push_back(tmp);
    }
    
    // Get the old step deltas for this function.  Convert it to a vector.
    const double* s=mxGetPr(prhs[1]);
    list < vector <double> > oldS;
    for(int j=0;j<noldS;j++){
    	vector <double> tmp(s+moldS*j,s+moldS*(j+1));
    	oldS.push_back(tmp);
    }

    // Extract the direction we evaluate the Hessian 
    const double* p0=mxGetPr(prhs[2]);
    vector <double> p(p0,p0+mp);
    
    // Allocate memory for the result of evaluation
    vector <double> result(mp);

    // Evaluate the BFGS preconditioner based on this data
    (Preconditioners::BFGS < vector <double> > (oldY,oldS))(p,result);

    // Create the output
    plhs[0] = mxCreateDoubleMatrix(mp,np,mxREAL);
    double* result0=mxGetPr(plhs[0]);
    for(int i=0;i<mp;i++)
    	result0[i]=result[i];
}
