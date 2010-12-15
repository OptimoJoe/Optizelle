#include "mex.h"
#include "peopt.h"
#include "simple_matching.h"

// Wrapper function for calling the example objective function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    // Check for valid inputs.  The order of inputs is
    // 1. A (cell array)
    // 2. B
    // 3. b (cell array)
    // 4. d
    // 5. u
    // 6. delta
    // 7. max_iter 
    // 8. eps_cg
    if(nrhs!=8)
	mexErrMsgTxt("The function getStep requires 8 inputs.");
    if(nlhs!=1)
	mexErrMsgTxt("The function getStep requires 1 output.");

    // Get the size of the collection of As
    const int mA=mxGetM(prhs[0]);
    const int nA=mxGetN(prhs[0]);

    // Get the size of B
    const int mB=mxGetM(prhs[1]);
    const int nB=mxGetN(prhs[1]);

    // Get the size of the collection of bs
    const int mb=mxGetM(prhs[2]);
    const int nb=mxGetN(prhs[2]);

    // Get the size of d
    const int md=mxGetM(prhs[3]);
    const int nd=mxGetN(prhs[3]);

    // Get the size of our control
    const int mu=mxGetM(prhs[4]);
    const int nu=mxGetN(prhs[4]);

    // Get the size of delta
    const int mdelta=mxGetM(prhs[5]);
    const int ndelta=mxGetN(prhs[5]);
    
    // Get the size of max_iter
    const int mmax_iter=mxGetM(prhs[6]);
    const int nmax_iter=mxGetN(prhs[6]);
    
    // Get the size of eps_cg 
    const double meps_cg=mxGetM(prhs[7]);
    const double neps_cg=mxGetN(prhs[7]);

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
    if(mB!=nB)
    	mexErrMsgTxt("The linear operators must be square.");
    if(mdelta!=1 && ndelta!=1)
    	mexErrMsgTxt("The parameter delta must be a scalar.");
    if(mmax_iter!=1 && nmax_iter!=1)
    	mexErrMsgTxt("The parameter max_iter must be a scalar.");
    if(meps_cg!=1 && neps_cg!=1)
    	mexErrMsgTxt("The parameter eps_cg must be a scalar.");

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

	// Verify this size is compatible with the measured data
	if(md!=mAAA)
	    mexErrMsgTxt("The measured data and the linear operators must have"
		" the same number of columns");

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

    // Get the data for this function.  Convert it to a vector.
    const double* d=mxGetPr(prhs[3]);
    vector < vector <double> > dd;
    for(int j=0;j<nd;j++){
    	vector <double> tmp(d+md*j,d+md*(j+1));
    	dd.push_back(tmp);
    }

    // Extract the control
    const double* uu=mxGetPr(prhs[4]);
    vector <double> u(uu,uu+nu);
    
    // Extract the initial trust region size 
    double delta=mxGetScalar(prhs[5]);
    
    // Extract the maximum number of iterations for Steihaug-Toint 
    const int max_iter=(int)mxGetScalar(prhs[6]);
    
    // Extract the stopping tolerance for truncated cg 
    const double eps_cg=mxGetScalar(prhs[7]);
    
    // Form the solution operator
    BasicOp h(A,B,b);

    // Create the matching operator
    BasicMatching f(dd);

    // Allocate memory for the workspaces
    list < vector <double> > workU(5);
    for(list < vector <double> >::iterator v=workU.begin();v!=workU.end();v++)
    	v->resize(nA);
    list < vector <double> > workU2(5);
    for(list < vector <double> >::iterator v=workU2.begin();v!=workU2.end();v++)
    	v->resize(nA);
    list < vector <double> > infoU(0);
    list < vector <double> > workY(5);
    for(list < vector <double> >::iterator v=workY.begin();v!=workY.end();v++)
    	v->resize(nB);

    // Allocate memory for the gradient
    vector <double> g(nA);

    // Allocate memory for the step
    vector <double> s(nA);

    // Find the gradient
    TrustRegion::getGradient < vector <double> , vector <double> >
	(f,h,workY,workU)(u,g);

    // Find the TR step using GaussNewton for the Hessian approximation with
    // no preconditioner
    typedef vector <double> U;
    typedef vector <double> Y;
    Preconditioners::Identity <U> Id;
    Hessians::GaussNewton <Y,U> GN (f,h,u,workY,workU2);
    TrustRegion::getStep <U> (Id,GN,u,g,delta,max_iter,eps_cg,workU,s);
    
    // Create the output
    plhs[0] = mxCreateDoubleMatrix(mu,nu,mxREAL);
    double* uuu=mxGetPr(plhs[0]);
    for(int i=0;i<nu;i++)
    	uuu[i]=s[i];
}
