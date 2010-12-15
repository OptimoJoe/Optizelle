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
    // 5. u0
    // 6. eta1
    // 7. eta2
    // 8. delta
    // 9. eps_cg
    // 10. eps_g
    // 11. eps_d
    // 12. max_iter for Steihaug-Toint
    // 13. H (Hessian) (0=Steepest Descent, 1=Gauss-Newton)
    // 14. Minv (Preconditioner) (0=Identity, 1=BFGS)
    if(nrhs!=14)
	mexErrMsgTxt("The function pe_test requires 14 inputs.");
    if(nlhs!=2)
	mexErrMsgTxt("The function pe_test requires 2 outputs.");

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

    // Get the size of eta1 
    const int meta1=mxGetM(prhs[5]);
    const int neta1=mxGetN(prhs[5]);

    // Get the size of eta2
    const int meta2=mxGetM(prhs[6]);
    const int neta2=mxGetN(prhs[6]);

    // Get the size of delta
    const int mdelta=mxGetM(prhs[7]);
    const int ndelta=mxGetN(prhs[7]);
    
    // Get the size of eps_cg 
    const int meps_cg=mxGetM(prhs[8]);
    const int neps_cg=mxGetN(prhs[8]);

    // Get the size of eps_g 
    const int meps_g=mxGetM(prhs[9]);
    const int neps_g=mxGetN(prhs[9]);

    // Get the size of eps_d 
    const int meps_d=mxGetM(prhs[10]);
    const int neps_d=mxGetN(prhs[10]);
    
    // Get the size of max_iter
    const int mmax_iter=mxGetM(prhs[11]);
    const int nmax_iter=mxGetN(prhs[11]);
    
    // Get the size of the Hessian approximation 
    const int mH=mxGetM(prhs[12]);
    const int nH=mxGetN(prhs[12]);
    
    // Get the size of the preconditioenr for truncated-CG 
    const int mMinv=mxGetM(prhs[13]);
    const int nMinv=mxGetN(prhs[13]);

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
    if(meta1!=1 && neta1!=1)
    	mexErrMsgTxt("The parameter eta1 must be a scalar.");
    if(meta2!=1 && neta2!=1)
    	mexErrMsgTxt("The parameter eta2 must be a scalar.");
    if(mdelta!=1 && ndelta!=1)
    	mexErrMsgTxt("The parameter delta must be a scalar.");
    if(meps_cg!=1 && neps_cg!=1)
    	mexErrMsgTxt("The parameter eps_cg must be a scalar.");
    if(meps_g!=1 && neps_g!=1)
    	mexErrMsgTxt("The parameter eps_g must be a scalar.");
    if(meps_d!=1 && neps_d!=1)
    	mexErrMsgTxt("The parameter eps_d must be a scalar.");
    if(mmax_iter!=1 && nmax_iter!=1)
    	mexErrMsgTxt("The parameter max_iter must be a scalar.");
    if(mH!=1 && nH!=1)
    	mexErrMsgTxt("The parameter H must be a scalar.");
    if(mMinv!=1 && nMinv!=1)
    	mexErrMsgTxt("The parameter Minv must be a scalar.");

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
    
    // Extract the first trust region parameter
    const double eta1=mxGetScalar(prhs[5]);
    
    // Extract the second trust region parameter
    const double eta2=mxGetScalar(prhs[6]);
    
    // Extract the initial trust region size 
    double delta=mxGetScalar(prhs[7]);

    // Extract the stopping tolerance for truncated cg 
    const double eps_cg=mxGetScalar(prhs[8]);

    // Extract the stopping tolerance for the gradient
    const double eps_g=mxGetScalar(prhs[9]);

    // Extract the stopping tolerance for the trust region size 
    const double eps_d=mxGetScalar(prhs[10]);
    
    // Extract the maximum number of iterations for Steihaug-Toint 
    const int max_iter=(int)mxGetScalar(prhs[11]);

    // Extract the Hessian approximation
    const int whichH=(int)mxGetScalar(prhs[12]);

    // Extract the preconditioner for truncated-CG 
    const int whichMinv=(int)mxGetScalar(prhs[13]);

    // Form the solution operator
    BasicOp h(A,B,b);

    // Create the matching operator
    BasicMatching f(dd);
   
    // Specify the domains for the control and state respectively
    typedef vector <double> U;
    typedef vector <double> Y;

    // Allocate memory for the workspaces
    list <Y> workHessianY;
    list <U> workHessianU;
    list <Y> workGradientY; 
    list <U> workGradientU; 
    list <Y> workObj;
    list <U> infoBFGS;
    list <U> workTR(4);
    for(list <U>::iterator v=workTR.begin(); v!=workTR.end(); v++)
	v->resize(nA);

    // Allocate memory for the gradient and the trial step
    vector <double> g(nA);
    vector <double> s(nA);

    // Allocate memory for the old gradient and trial step
    vector <double> g_old(nA);
    vector <double> s_old(nA);

    // Set which Hessian approximation we'll use
    typedef vector <double> Y;
    typedef vector <double> U;
    Operator <U,U>* H;
    switch(whichH){
   	case 0:
	    H=new Hessians::Identity<U>();
	    break;
	case 1:
	    // In the case of Gauss-Newton, we require three work elements
	    // in the state space and 1 in the control
	    workHessianU.resize(1);
	    for(list <U>::iterator v=workHessianU.begin();
		v!=workHessianU.end();
		v++
	    ) v->resize(nA);

	    workHessianY.resize(3);
	    for(list <Y>::iterator v=workHessianY.begin();
		v!=workHessianY.end();
		v++
	    ) v->resize(nB);

	    // Create the approximation. NOTE: this is parameterized on u
	    H=new Hessians::GaussNewton<Y,U>(f,h,u,workHessianY,workHessianU);
	    break;
	case 2:
	    // In the case of full-Newton, we require four work elements
	    // in the state space and 1 in the control 
	    workHessianU.resize(1);
	    for(list <U>::iterator v=workHessianU.begin();
		v!=workHessianU.end();
		v++
	    ) v->resize(nA);

	    workHessianY.resize(4);
	    for(list <Y>::iterator v=workHessianY.begin();
		v!=workHessianY.end();
		v++
	    ) v->resize(nB);

	    // Create the Hessian approximation
	    H=new Hessians::Newton<Y,U>(f,h,u,workHessianY,workHessianU);
	    break;
	default:
	    mexErrMsgTxt("Invalid Hessian selection.");
	    break;
    }

    // Set which preconditioner we'll use
    Operator <U,U>* Minv;
    switch(whichMinv){
    	case 0:
	    Minv=new Preconditioners::Identity <U>();
	    break;
	case 1:
	    Minv=new Preconditioners::BFGS <U>(infoBFGS);
	    break;
	default:
	    mexErrMsgTxt("Invalid preconditioner selection.");
	    break;	
    }

    // Find the objective value operator
    workObj.resize(2);
    for(list <Y>::iterator v=workHessianY.begin(); v!=workHessianY.end(); v++)
	v->resize(nB);
    
    TrustRegion::getObjValue <Y,U> getObjValue(f,h,workObj);
    	
    // Find the gradient operator
    workGradientU.resize(1);
    for(list <U>::iterator v=workGradientU.begin();
	v!=workGradientU.end();
	v++
    ) v->resize(nA);

    workGradientY.resize(3);
    for(list <Y>::iterator v=workGradientY.begin();
	v!=workGradientY.end();
	v++
    ) v->resize(nB);

    TrustRegion::getGradient <Y,U>
	getGradient(f,h,workGradientY,workGradientU);
   

    // Find the gradient
    getGradient(u,g);

    // Track the number of iterations required to converge
    int iters=0;

    // While the stopping criteria is not satisfied
    while(!TrustRegion::checkStop <U> (*Minv,u,g,delta,eps_g,eps_d,workTR)){

	// Find the new step
	TrustRegion::getStep <U> (*Minv,*H,u,g,delta,max_iter,eps_cg,workTR,s);

	// Check whether this step provides sufficient reduction
	bool accept;
	double obj_u;
	double obj_ups;
	double rho;
	TrustRegion::checkStep <U> (u,s,getObjValue,*H,g,
	    eta1,eta2,workTR,delta,accept,obj_u,obj_ups,rho);

	// If we're not on the first iteration and we accept the step,
	// store the old gradient and trial step.
	if(iters!=0 && accept){
	    Operations::copy(g,g_old);
	    Operations::copy(s,s_old);
	}

	// If the step is acceptable, move to the new point
	if(accept) Operations::axpy(1.0,s,u);
    	
	// Find the gradient
	getGradient(u,g);
	
	// If we're not on the first iteration and we accept the step, find
	// y_k = grad f(h(u_k)) - grad f(h(u_{k-1}))
	// s_k = u_k-u_{k-1}
	if(iters!=0 && accept && whichMinv==1){
	    vector <double>& tmp=*(workTR.begin());
	    Operations::copy(s,tmp);
	    Operations::axpy(-1.,s_old,tmp);
	    infoBFGS.push_front(tmp);

	    Operations::copy(g,tmp);
	    Operations::axpy(-1.,g_old,tmp);
	    infoBFGS.push_front(tmp);
	}

	// Increment the iteration count
	iters++;
    }

    // Create the output
    // The control resulting from inversion
    plhs[0] = mxCreateDoubleMatrix(mu,nu,mxREAL);
    double* uuu=mxGetPr(plhs[0]);
    for(int i=0;i<nu;i++)
    	uuu[i]=u[i];

    // The number of iterations required to converge
    const int dims[1]={1};
    plhs[1] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    int* niters=(int*)mxGetData(plhs[1]);
    niters[0]=iters;
}
