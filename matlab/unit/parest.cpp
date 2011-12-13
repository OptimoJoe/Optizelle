#include "mex.h"
#include "core.h"
#include "simple_matching.h"

#include <vector>
#include <string>
#include <iostream>
using namespace std;

void pe_error(const char* msg){
    mexErrMsgTxt(msg);
}

void optimize(
    vector < vector<double> >& A,
    vector <double>& B,
    vector < vector<double> >& b,
    vector < vector<double> >& d,
    vector <double>& u,
    string& H_type,
    string& Minv_type,
    double eps_g,
    double eps_d,
    double eps_f,
    double eps_cg,
    int max_cg_iter,
    double eta1,
    double eta2,
    double delta,
    int stored_history,
    int history_reset,
    int max_iter
){
    // Create a shortcut for the control and state space
    typedef vector <double> U;
    typedef vector <double> Y;

    // Set the precision for cout
    cout.precision(8);
    scientific(cout);

    // Create the matching operator
    BasicMatching f(d);

    // Create the solution operator
    BasicOp h(A,B,b); 

    // Create the objective function operator
    list <Y> workObj;
    for(int i=0;i<2;i++){
    	Y tmp(d[0]);
	workObj.push_back(tmp);
    }
    TrustRegion::getObjValue <Y,U> ObjFunc(f,h,workObj);

    // Create the gradient operator
    list <Y> workGrad_Y;
    for(int i=0;i<3;i++){
    	Y tmp(d[0]);
	workGrad_Y.push_back(tmp);
    }
    list <U> workGrad_U;
    for(int i=0;i<1;i++){
    	Y tmp(u);
	workGrad_U.push_back(tmp);
    }
    TrustRegion::getGradient <Y,U> Gradient(f,h,workGrad_Y,workGrad_U);

    // Allocate memory for the gradient
    U g(u);

    // Allocate memory for a typical gradient
    U g_typ(u);

    // Allocate memory for a precondition gradient we use in truncated cg
    U g_tcg(u);

    // Allocate memory for the trial step
    U s(u);

    // Allocate memory so that we can store an old gradient, step,
    // and a list of vales used in quasi-Newton methods
    U g_old(u);
    U u_old(u);
    list <U> oldY;
    list <U> oldS;

    // Allocate some extra memory for the Hessian approximation if needed
    list <U> workH_U;
    list <U> workH_Y;

    // Allocate some extra memory for the preconditioner if needed
    list <U> workMinv;

    // Allocate memory for the workspaces
    list <U> workU(5);
    for(list <U>::iterator v=workU.begin();v!=workU.end();v++)
        v->resize(u.size());
    list < vector <double> > workY(5);
    for(list <Y>::iterator v=workY.begin();v!=workY.end();v++)
        v->resize(b[0].size());

    // Create the Hessian approximation
    Operator <U,U>* H;
    if(H_type=="Identity")
    	H=new Hessians::Identity<U>();
    else if(H_type=="BFGS")
    	H=new Hessians::BFGS<U>(oldY,oldS,workH_U);
    else if(H_type=="SR1")
    	H=new Hessians::SR1<U>(oldY,oldS,workH_U);
    else if(H_type=="GaussNewton"){
    	for(int i=0;i<3;i++){
	    Y tmp(d[0]);
	    workH_Y.push_back(tmp);
	}
	U tmp(u);
	workH_U.push_back(tmp);
    	H=new Hessians::GaussNewton <Y,U> (f,h,u,workH_Y,workH_U);
    } else if(H_type=="Newton"){
    	for(int i=0;i<4;i++){
	    Y tmp(d[0]);
	    workH_Y.push_back(tmp);
	}
	U tmp(u);
	workH_U.push_back(tmp);
    	H=new Hessians::Newton <Y,U> (f,h,u,workH_Y,workH_U);
    } else 
    	pe_error("Detected an unknown Hessian approximation.");

    // Create the preconditioner
    Operator <U,U>* Minv;
    if(Minv_type=="Identity")
    	Minv=new Preconditioners::Identity<U>();
    else if(Minv_type=="BFGS")
    	Minv=new Preconditioners::BFGS<U>(oldY,oldS);
    else if(Minv_type=="SR1")
    	Minv=new Preconditioners::SR1<U>(oldY,oldS,workMinv);
    else 
    	pe_error("Detected an unknown preconditioner.");

    // Check if we have, in theory, a perfect preconditioner
    bool perfect_preconditioner;
    if( (Minv_type=="BFGS" && H_type=="BFGS") ||
	(Minv_type=="SR1" && H_type=="SR1"))
    	perfect_preconditioner=true;
    else
    	perfect_preconditioner=false;

    // Keep track of the number of problem iterations
    int iter=1;

    // Track the number of times we reject a step and the trust region
    // radius prior to the search for a trial step
    int nreject=0;
    double delta0=delta;

    // Evaluate the gradient and objective function prior to beginning the
    // optimization
    cout << "Computing gradient before optimization begins." << endl;
    Gradient(u,g);
    if(perfect_preconditioner) (*Minv)(g,g_tcg);

    cout << "Computing the objective function before optimization begins."
	<< endl;
    double obj_u=ObjFunc(u);
    double obj_typ=obj_u;
    double obj_ums=std::numeric_limits <double> ::quiet_NaN();

    cout << "Initial objective value:\t\t\t" << obj_u << endl;
    double inorm_g = sqrt(Operations::innr(g,g));
    cout << "Initial norm of the gradient:\t\t\t" << inorm_g <<endl;

    // While the stopping criteria is not satisfied
    TrustRegion::StoppingCondition reason;

    while((reason=TrustRegion::checkStop <U> (*Minv,u,g,g_typ,obj_u,obj_ums,
	delta,eps_g,eps_d,eps_f,iter,max_iter,workU))==TrustRegion::NotConverged
    ){
	// Print out some iteration information
	cout << endl << "Optimization iteration:\t\t\t" << iter << endl;
	  
	// Output the current trust region radius
	cout << "Current trust region radius:\t\t\t" << delta << endl;

	// Find the new step
	cout << "Running truncated CG to find a new trial step." << endl;
       
	// If we have a perfect preconditioner, simplify the call and use
	// the preconditioned gradient
	if(perfect_preconditioner){
	    Preconditioners::Identity<U> P_Id;
	    Hessians::Identity<U> H_Id;
	    TrustRegion::getStep <U>
	      (P_Id,H_Id,u,g_tcg,delta,max_cg_iter,eps_cg,workU,s);

	// Otherwise, calculate the step as normal
	} else
	    TrustRegion::getStep <U>
		(*Minv,*H,u,g,delta,max_cg_iter,eps_cg,workU,s);

	// Check whether this step provides sufficient reduction
	bool accept;
	double obj_ups;
	double rho;
	cout << "Calculating reduction in the objective." << endl;
	TrustRegion::checkStep <U> (u,s,ObjFunc,*H,g,obj_u,eta1,eta2,workU,
	  delta,accept,obj_ups,rho);

	// Output some diagnostic information about the TR process.  Both
	// of these numbers can be NaNs, so in order to fix the
	// formatting issues we check whether or not the number equals
	// itself which is false for NaN
	cout << "Objective at current step:\t\t\t" << obj_u << endl;
	if(obj_ups==obj_ups)
	    cout << "Objective at trial step:\t\t\t" << obj_ups << endl;
	else
	    cout << "Objective at trial step:\t\t\t" << obj_ups << endl;
	if(rho==rho)
	    cout << "Ratio between predicted and actual reduction:\t" << rho
		<< endl;
	else
	    cout << "Ratio between predicted and actual reduction:\t" << rho
		<< endl;

#if 0
	// If we accept a step, do a finite difference check on the direction
	if(accept){
	  cout << "Directional derivative in the direction of the accepted "
	    "step." << endl;
	  
	  // Begin begin calculating the directional derivative via the
	  // gradient
	  double dd_grad=Operations::innr(g,s);

	  // Next, do a two point finite difference using a variety of
	  // step lengths
	  U& work1=*(workU.begin());
	  for(int i=1;i<=10;i++){
	    double epsilon=pow(.1,i);
	    Operations::copy(u,work1);
	    Operations::axpy(epsilon,s,work1);
	    double obj_upes=ObjFunc(work1);
	    double dd=(obj_upes-obj_u)/epsilon;
	    cout << "The relative difference (1e-" << i << "): "
	      << fabs(dd_grad-dd)/(1+fabs(dd_grad)) << endl;
	  }
	}
#endif

	// Actions to take if we accept a step
	if(accept){
	    // Output some diagnostic information bout this objective value
	    cout << "Accepted trial step:\t\t\t\ttrue" << endl; 

	    // Before we move to the new trial step and calculate the new
	    // gradient, store the old step and old gradient for quasi-Newton
	    // methods
	    Operations::copy(u,u_old);
	    Operations::copy(g,g_old);

	    // Actually move to the new trial step
	    Operations::axpy(1.0,s,u);
	    
	    // Compute the gradient
	    cout << "Computing the gradient." << endl;
	    Gradient(u,g);
	    
	    // Next, store the objective value at the new trial step so
	    // we don't have to recompute it
	    obj_ums=obj_u;
	    obj_u=obj_ups;
	    cout << "Current objective value:\t\t\t" << obj_u << endl;
	    
	    // If we are storing our history and we are not on the first
	    // iteration, store the difference between 
	    // the previous step and this one as well as the previous gradient
	    // and this one.  In other words,
	    // y_k = grad f(h(u_k)) - grad f(h(u_{k-1}))
	    // s_k = u_k-u_{k-1}
	    if(stored_history>0){
	      // Determine s_k
	      U& s_k=*(workU.begin());
	      Operations::copy(u,s_k);
	      Operations::axpy(-1.,u_old,s_k);

	      // Determine y_k
	      U& y_k=*(++workU.begin());
	      Operations::copy(g,y_k);
	      Operations::axpy(-1.,g_old,y_k);

	      // If we're doing BFGS, make sure that <y_k,s_k> > 0.  Otherwise,
	      // don't store the vector
	      double inner_yk_sk=Operations::innr(y_k,s_k);
	      if(H_type!="BFGS" || inner_yk_sk > 0){
		  // If we're at the limit of the number of stored vectors,
		  // eliminate the last stored gradient and point
		  if(oldY.size()==stored_history){
		    oldY.pop_back();
		    oldS.pop_back();
		  }

		  // Store the resulting s_k
		  oldS.push_front(s_k);

		  // Store the resulting y_k
		  oldY.push_front(y_k);

		  // Determine if we need more work space
		  if(H_type=="BFGS" || H_type=="SR1"){
		      if(workH_U.size()<stored_history){
			  U tmp(u);
			  workH_U.push_front(tmp);
		      }
		  }
		  if(Minv_type=="SR1"){
		      if(workMinv.size()<stored_history){
			  U tmp(u);
			  workMinv.push_front(tmp);
		      }
		  }
	      }
	  }
	  
	  // After we update the quasi-Newton method, precondition the
	  // gradient if necessary.
	  if(perfect_preconditioner) (*Minv)(g,g_tcg);

	  // Output the norm of the gradient
	  double norm_g = sqrt(Operations::innr(g,g));
	  cout << "Norm of the gradient:\t\t\t\t" << norm_g << endl;

	  // Increment the iteration count
	  iter++;

	  // Reset the rejection counter and stored trust region radius
	  nreject=0;
	  delta0=delta;

	// If we reject the step
	} else{
	  // Keep track of how many times we reject a step
	  nreject++;

	  // If we reject a step too many times, reset the quasi-Newton
	  // approximation
	  if(history_reset!=0 && nreject>=history_reset && oldY.size() > 0){
	    // Eliminate all the stored history
	    oldY.resize(0);
	    oldS.resize(0);
	    if( H_type=="BFGS" || H_type=="SR1")
	      workH_U.resize(0);
	    if(Minv_type=="SR1")
	      workMinv.resize(0);

	    // Reset the trust region radius
	    delta=delta0;

	    // If we're using perfect preconditioners, make sure that we
	    // recalculate the gradient
	    if(perfect_preconditioner) Operations::copy(g,g_tcg);

	    cout << "Resetting the quasi-Newton approximation." << endl;
	  }
	  
	  // Output some diagnostic information about the objective value
	  cout << "Accepted trial step:\t\t\t\tfalse" << endl << endl; 
	}
    }

    // Output some diagnostic information about the current solution
    // and why we have converged 
    cout << endl << "The method has converged." << endl;
    cout << "The reason for convergence is:\t";
    switch(reason){
    case TrustRegion::GradientSmall:
      cout << "Norm of the gradient < eps_g" << endl;
      break;
    case TrustRegion::RelativeGradientSmall:
      cout << "Relative norm of the gradient < eps_g" << endl;
      break;
    case TrustRegion::TrustRegionSmall:
      cout << "Size of the trust region radius < eps_d" << endl;
      break;
    case TrustRegion::RelativeTrustRegionSmall:
      cout<<"Relative size of the trust region radius < eps_d"<<endl;
      break;
    case TrustRegion::RelativeObjectiveSmall:
      cout << "Relative change in objective value < eps_f" << endl;
      break;
    case TrustRegion::MaxItersExceeded:
      cout << "Number of optimization iterations > max_iter" << endl;
      break;
    case TrustRegion::NotConverged:
      cout << "Number of subiterations > max_sub_iter" << endl;
      break;
    }

    double norm_g = sqrt(Operations::innr(g,g));
    double norm_gtyp = sqrt(Operations::innr(g_typ,g_typ));
    cout << "Initial norm of the gradient:\t\t" << norm_gtyp << endl;
    cout << "Current norm of the gradient:\t\t" << norm_g << endl;
    cout << "Initial objective value:\t\t" << obj_typ << endl;
    cout << "Current objective value:\t\t" << obj_u << endl;
    cout << "Current trust region radius:\t\t" << delta << endl;
    cout << "Relative change in objective:\t\t" 
      << (obj_ums-obj_u)/(1.+fabs(obj_ums)) << endl << endl;
}

  // Wrapper function for calling the optimization 
  void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  // Check for valid inputs.  The order of inputs is
  // 1. A (cell array)
  // 2. B
  // 3. b (cell array)
  // 4. d (cell array)
  // 5. u0
  // 6. Hessian approximation (string)
  // 7. Preconditioner (string)
  // 8. eps_g
  // 9. eps_d
  // 10. eps_f
  // 11. eps_cg
  // 12. max_cg_iter
  // 13. eta1
  // 14. eta2
  // 15. delta
  // 16. stored_history
  // 17. history_reset
  // 18. max_iter
  if(nrhs!=18)
      mexErrMsgTxt("The function parest requires 18 inputs.");
  if(nlhs!=1)
      mexErrMsgTxt("The function parest requires 1 output.");

  // Get the size of the collection of As
  const int mA=mxGetM(prhs[0]);
  const int nA=mxGetN(prhs[0]);

  // Get the size of B
  const int mB=mxGetM(prhs[1]);
  const int nB=mxGetN(prhs[1]);

  // Get the size of the collection of bs
  const int mb=mxGetM(prhs[2]);
  const int nb=mxGetN(prhs[2]);

  // Get the size of the collection of ds
  const int md=mxGetM(prhs[3]);
  const int nd=mxGetN(prhs[3]);

  // Get the size of our initial control
  const int mu=mxGetM(prhs[4]);
  const int nu=mxGetN(prhs[4]);

  // Get the size of the Hessian approximation
  const int mH=mxGetM(prhs[5]);
  const int nH=mxGetN(prhs[5]);

  // Get the size of the preconditioner
  const int mMinv=mxGetM(prhs[6]);
  const int nMinv=mxGetN(prhs[6]);

  // Get the size of the gradient tolerance
  const int meps_g=mxGetM(prhs[7]);
  const int neps_g=mxGetN(prhs[7]);

  // Get the size of the trust region tolerance
  const int meps_d=mxGetM(prhs[8]);
  const int neps_d=mxGetN(prhs[8]);

  // Get the size of the objective change tolerance
  const int meps_f=mxGetM(prhs[9]);
  const int neps_f=mxGetN(prhs[9]);

  // Get the size of the truncated cg tolerance
  const int meps_cg=mxGetM(prhs[10]);
  const int neps_cg=mxGetN(prhs[10]);

  // Get the size of the maximum number of CG iterations 
  const int mmax_cg_iter=mxGetM(prhs[11]);
  const int nmax_cg_iter=mxGetN(prhs[11]);

  // Get the size of the lower TR tolerance 
  const int meta1=mxGetM(prhs[12]);
  const int neta1=mxGetN(prhs[12]);

  // Get the size of the upper TR tolerance 
  const int meta2=mxGetM(prhs[13]);
  const int neta2=mxGetN(prhs[13]);

  // Get the size of the initial trust region
  const int mdelta=mxGetM(prhs[14]);
  const int ndelta=mxGetN(prhs[14]);

  // Get the size of the stored history
  const int mstored_history=mxGetM(prhs[15]);
  const int nstored_history=mxGetN(prhs[15]);

  // Get the size of the history reset
  const int mhistory_reset=mxGetM(prhs[16]);
  const int nhistory_reset=mxGetN(prhs[16]);

  // Get the size of the maximum number of iterations
  const int mmax_iter=mxGetM(prhs[17]);
  const int nmax_iter=mxGetN(prhs[17]);

  // Verify the sizes
  if(mA!=1)
      mexErrMsgTxt("The cell array containing the linear operators must "
	  "contain only a single row.");
  if(mB!=nB)
      mexErrMsgTxt("The constant linear operator must be square.");
  if(mb!=1)
      mexErrMsgTxt("The cell array containing the sources must contain only "
	  "a single row.");
  if(md!=1)
      mexErrMsgTxt("The cell array containing the data must contain only "
	  "a single row.");
  if(nu!=1)
      mexErrMsgTxt("The initial control must contain a single column.");
  if(mu!=nA)
      mexErrMsgTxt("The size of the initial control must match the number of "
	  "linear operators.");
  if(mxIsChar(prhs[5])!=1)
      mexErrMsgTxt("The Hessian approximation must be a string.");
  if(mH!=1)
      mexErrMsgTxt("The Hessian approximation must be a string, row vector.");
  if(mxIsChar(prhs[6])!=1)
      mexErrMsgTxt("The preconditioner must be a string.");
  if(mMinv!=1)
      mexErrMsgTxt("The preconditioner must be a string, row vector.");
  if(meps_g!=1 && neps_g!=1)
      mexErrMsgTxt("The tolerance eps_g must be a scalar.");
  if(meps_d!=1 && neps_d!=1)
      mexErrMsgTxt("The tolerance eps_d must be a scalar.");
  if(meps_f!=1 && neps_f!=1)
      mexErrMsgTxt("The tolerance eps_f must be a scalar.");
  if(meps_cg!=1 && neps_cg!=1)
      mexErrMsgTxt("The tolerance eps_cg must be a scalar.");
  if(mmax_cg_iter!=1 && nmax_cg_iter!=1)
      mexErrMsgTxt("The tolerance max_cg_iter must be a scalar.");
  if(meta1!=1 && neta1!=1)
      mexErrMsgTxt("The tolerance eta1 must be a scalar.");
  if(meta2!=1 && neta2!=1)
      mexErrMsgTxt("The tolerance eta2 must be a scalar.");
  if(mdelta!=1 && ndelta!=1)
      mexErrMsgTxt("The initial trust region radius must be a scalar.");
  if(mstored_history!=1 && nstored_history!=1)
      mexErrMsgTxt("The amount of stored history must be a scalar.");
  if(mhistory_reset!=1 && nhistory_reset!=1)
      mexErrMsgTxt("The frequency that we reset the history must be "
	  "a scalar.");
  if(mmax_iter!=1 && nmax_iter!=1)
      mexErrMsgTxt("The maximum number of iterations must be a scalar.");

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

  // Extract the linear operator that doesn't depend on the control
  const double* BB=mxGetPr(prhs[1]);
  vector <double> B(BB,BB+mB*nB);

  // Extract the sources 
  if(!mxIsCell(prhs[2]))
      mexErrMsgTxt("The third argument must be a cell array containing the "
	  "sources for the problem.");
  vector < vector<double> > b;
  for(int i=0;i<nb;i++){
      // Extract the cell
      mxArray* bbb=mxGetCell(prhs[2],i);

      // Find its size
      const int mbbb=mxGetM(bbb);
      const int nbbb=mxGetN(bbb);

      // Verify this size aligns with B and has one column
      if(mbbb!=mB || nbbb!=1)
	  mexErrMsgTxt("The size of the sources must align with B and have "
	      "a single column.");

      // Extract a double pointer from this matrix
      const double* bb=mxGetPr(bbb);
      vector <double> tmp(bb,bb+mbbb);
      b.push_back (tmp);
  }

  // Extract the data 
  if(!mxIsCell(prhs[3]))
      mexErrMsgTxt("The fourth argument must be a cell array containing the "
	  "sources for the problem.");
  vector < vector<double> > d;
  for(int i=0;i<nd;i++){
      // Extract the cell
      mxArray* ddd=mxGetCell(prhs[3],i);

      // Find its size
      const int mddd=mxGetM(ddd);
      const int nddd=mxGetN(ddd);

      // Verify this size aligns with B and has one column
      if(mddd!=mB || nddd!=1)
	  mexErrMsgTxt("The size of the data must align with B and have "
	      "a single column.");

      // Extract a double pointer from this matrix
      const double* dd=mxGetPr(ddd);
      vector <double> tmp(dd,dd+mddd);
      d.push_back (tmp);
  }

  // Extract the control
  const double* uu=mxGetPr(prhs[4]);
  vector <double> u(uu,uu+mu);

  // Extract the Hessian approximation
  const char* HH;
  HH=mxArrayToString(prhs[5]);
  string H(HH);

  // Extract the preconditioner
  const char* MMinv;
  MMinv=mxArrayToString(prhs[6]);
  string Minv(MMinv);

  // Extract several tolerances and parameters
  double eps_g=mxGetScalar(prhs[7]);
  double eps_d=mxGetScalar(prhs[8]);
  double eps_f=mxGetScalar(prhs[9]);
  double eps_cg=mxGetScalar(prhs[10]);
  int max_cg_iter=(int)mxGetScalar(prhs[11]);
  double eta1=mxGetScalar(prhs[12]);
  double eta2=mxGetScalar(prhs[13]);
  double delta=mxGetScalar(prhs[14]);
  int stored_history=(int)mxGetScalar(prhs[15]);
  int history_reset=(int)mxGetScalar(prhs[16]);
  int max_iter=(int)mxGetScalar(prhs[17]);

  // Do some more input validation
  if(eps_g<=0)
      mexErrMsgTxt("The parameter eps_g must a be positive real.");
  if(eps_d<=0)
      mexErrMsgTxt("The parameter eps_d must a be positive real.");
  if(eps_f<=0)
      mexErrMsgTxt("The parameter eps_f must a be positive real.");
  if(eps_cg<=0)
      mexErrMsgTxt("The parameter eps_cg must a be positive real.");
  if(max_cg_iter<=0)
      mexErrMsgTxt("The parameter max_cg_iter must be positive integer.");
  if(eta1<=0)
      mexErrMsgTxt("The parameter eta1 must a be positive real.");
  if(eta2<=0)
      mexErrMsgTxt("The parameter eta2 must a be positive real.");
  if(delta<=0)
      mexErrMsgTxt("The initial trust region radius must "
	  "be a positive real.");
  if(stored_history<0)
      mexErrMsgTxt("The amount of stored history must be a "
	  "nonnegative integer.");
  if(history_reset<0)
      mexErrMsgTxt("The safety tolerance for resetting the history must "
	  "be a nonnegative integer");
  if(max_iter<0)
      mexErrMsgTxt("The maximum number of iterations must "  
	  "be a nonnegative integer");

  // Run the optimization
  optimize(A,B,b,d,u,H,Minv,eps_g,eps_d,eps_f,eps_cg,max_cg_iter,eta1,eta2,
      delta,stored_history,history_reset,max_iter);

  // Create the output
  plhs[0] = mxCreateDoubleMatrix(mu,1,mxREAL);
  double* uu_out=mxGetPr(plhs[0]);
  for(int i=0;i<mu;i++)
      uu_out[i]=u[i];
  }

