#ifndef PELAB_H 
#define PELAB_H 

#include "mex.h"
#include <string>
#include <iostream>
#include <memory>
#include "peopt.h"

struct pelab{
    struct VS{
	typedef double Real;
	typedef mxArray* Vector;

	static const std::string copy_,scal_,axpy_,innr_,init_;
	
	// y <- x (Shallow.  No memory allocation.)
	static void copy(const Vector& x, Vector& y){
	    // Create the inputs and outputs
	    Vector input[1]={x};

	    // Call the copy function from matlab
	    mexCallMATLAB(1,&y,1,input,copy_.c_str());
	}

	// x <- alpha * x
	static void scal(const Real& alpha_, Vector& x){
	    // Create memory for the scalar
	    mxArray* alpha;
	    alpha=mxCreateDoubleMatrix(1,1,mxREAL);
	    mxGetPr(alpha)[0]=alpha_;

	    // Create the inputs and outputs 
	    Vector input[2]={alpha,x};

	    // Compute the scalar multiplication
	    mexCallMATLAB(1,&x,2,input,scal_.c_str());

	    // Free memory from the scalar
	    mxDestroyArray(alpha);
	}

	// y <- alpha * x + y
	static void axpy(const Real& alpha_, const Vector& x, Vector& y){
	    // Create memory for the scalar
	    mxArray* alpha;
	    alpha=mxCreateDoubleMatrix(1,1,mxREAL);
	    mxGetPr(alpha)[0]=alpha_;

	    // Create the inputs and outputs 
	    Vector input[3]={alpha,x,y};

	    // Compute the addition 
	    mexCallMATLAB(1,&y,3,input,axpy_.c_str());

	    // Free memory from the scalar
	    mxDestroyArray(alpha);
	    
	}

	// innr <- <x,y>
	static Real innr(const Vector& x,const Vector& y){
	    // Create memory for the result 
	    mxArray* alpha;
	    alpha=mxCreateDoubleMatrix(1,1,mxREAL);
	    
	    // Create the inputs and outputs 
	    Vector input[2]={x,y};
	    
	    // Compute the inner product 
	    mexCallMATLAB(1,&alpha,2,input,innr_.c_str());

	    // Get the result of the computation
	    double alpha_=mxGetScalar(alpha);

	    // Free memory from the scalar
	    mxDestroyArray(alpha);

	    // Return the result
	    return alpha_;
	}

	// Memory allocation and size setting
	static void init(const Vector& x, Vector& y){
	    // Create the inputs and outputs
	    Vector input[1]={x};

	    // Call the init function from matlab
	    mexCallMATLAB(1,&y,1,input,init_.c_str());
	}

	// Error message
	static void error(std::string msg){
	    mexErrMsgTxt(msg.c_str());
	}

	// Print message
	static void print(std::string msg){
	    std::cout << msg;	
	}
    };

    typedef peopt::DataStructures <VS>::List VectorList;
	
    class Functional: public peopt::Functional <VS> {
    private:
    	mxArray*& f;
    public:
    	Functional(mxArray*& f_) : f(f_) {}

	// Basic application
	VS::Real operator () (const VS::Vector& x) const {
	    // Create memory for the result 
	    mxArray* alpha;
	    alpha=mxCreateDoubleMatrix(1,1,mxREAL);
	    
	    // Create the inputs 
	    mxArray* input[2]={f,x};
	    
	    // Compute the objective function
	    mexCallMATLAB(1,&alpha,2,input,"feval");

	    // Get the result of the computation
	    double alpha_=mxGetScalar(alpha);

	    // Free memory from the scalar
	    mxDestroyArray(alpha);

	    // Return the result
	    return alpha_;
	}
    };
    
    class Gradient: public peopt::Operator <VS,VS> {
    private:
    	mxArray*& g;
    public:
    	Gradient(mxArray*& g_) : g(g_) {}

	// Basic application
	void operator () (const VS::Vector& x,VS::Vector& y) const {
	    // Create the inputs 
	    mxArray* input[2]={g,x};
	    
	    // Compute the gradient function
	    mexCallMATLAB(1,&y,2,input,"feval");
	}
    };
    
    class Operator: public peopt::Operator <VS,VS> {
    private:
    	mxArray*& op;
	VS::Vector& base;  // The point around which we find the operator 
    public:
    	Operator(mxArray*& op_,VS::Vector& base_)
	    : op(op_), base(base_) {}

	// Basic application
	void operator () (const VS::Vector& x,VS::Vector& y) const {
	    // Create the inputs 
	    mxArray* input[3]={op,base,x};
	    
	    // Compute the gradient function
	    mexCallMATLAB(1,&y,3,input,"feval");
	}
    };

    class StateManipulator : public peopt::core <VS>::StateManipulator {
    private:
    	mxArray*& f;
    public:
    	StateManipulator(mxArray*& f_) : f(f_) {}

    	void operator () (peopt::core <VS>::State& state){
	    // Create structs for the original and manipulated state 
	    mxArray* original;
	    mxArray* manipulated;

	    // Copy the state to the matlab structure
	    stateToStruct(state,original);

	    // Allocate memory for the manipulated state
	    manipulated=mxDuplicateArray(original);

	    // Setup the inputs for the manipulator
	    mxArray* input[2]={f,original};

	    // Allow the user to manipulate the state 
	    mexCallMATLAB(1,&manipulated,2,input,"feval");
	    //mexCallMATLAB(0,NULL,1,&original,fname.c_str());

	    // Copy in the user changes
	    structToState(manipulated,state);

	    // Free the matlab memory 
	    //mxFree(original);
	    //mxFree(manipulated);
	}
    };


    // Converts a peopt state to a structure usable by matlab
    static void stateToStruct(
	peopt::core <VS>::State& state,
	mxArray*& mstate,
	const bool brief=false
    ){ 
	// Create a structure with all the parameters
	const char* params[44]={
	    "eps_g",
	    "eps_s",
	    "stored_history",
	    "history_reset",
	    "iter",
	    "iter_max",
	    "opt_stop",
	    "krylov_iter",
	    "krylov_iter_max",
	    "krylov_iter_total",
	    "krylov_stop",
	    "krylov_rel_err",
	    "eps_krylov",
	    "algorithm_class",
	    "Minv_type",
	    "H_type",
	    "norm_g",
	    "norm_gtyp",
	    "norm_s",
	    "norm_styp",
	    "u",
	    "g",
	    "s",
	    "u_old",
	    "g_old",
	    "s_old",
	    "oldY",
	    "oldS",
	    "obj_u",
	    "obj_ups",
	    "verbose",
	    "delta",
	    "delta_max",
	    "eta1",
	    "eta2",
	    "rho",
	    "rejected_trustregion",
	    "alpha",
	    "linesearch_iter",
	    "linesearch_iter_max",
	    "linesearch_iter_total",
	    "eps_ls",
	    "dir",
	    "kind"};

	const char* params_brief[25]={
            "eps_g",
	    "eps_s",
	    "stored_history",
	    "history_reset",
            "iter_max",
	    "krylov_iter_max",
	    "eps_krylov",
	    "algorithm_class",
            "Minv_type",
	    "H_type",
	    "verbose",
	    "delta",
	    "delta_max",
	    "eta1",
	    "eta2",
            "alpha",
	    "linesearch_iter_max",
	    "eps_ls",
	    "dir",
	    "kind",
            "F","G","H","Minv","StateManipulator"};

	if(!brief)
	    mstate=mxCreateStructMatrix(1,1,44,params);
	else
	    mstate=mxCreateStructMatrix(1,1,25,params_brief);

	// Fill the parameters based on what is found in the state
	mxSetField(mstate,0,"eps_g",mxCreateDoubleScalar(state.eps_g)); 
	mxSetField(mstate,0,"eps_s",mxCreateDoubleScalar(state.eps_s)); 

	mxSetField(mstate,0,"stored_history",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"stored_history")))[0]
	    =state.stored_history;

	mxSetField(mstate,0,"history_reset",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"history_reset")))[0]
	    =state.history_reset;

	mxSetField(mstate,0,"iter_max",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"iter_max")))[0]=state.iter_max;
	
	mxSetField(mstate,0,"krylov_iter_max",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"krylov_iter_max")))[0]
	    =state.krylov_iter_max;

	mxSetField(mstate,0,"eps_krylov",
	    mxCreateDoubleScalar(state.eps_krylov)); 

	switch(state.algorithm_class){
	    case peopt::TrustRegion:
		mxSetField(mstate,0,"algorithm_class",
		    mxCreateString("TrustRegion"));
		break;
	    case peopt::LineSearch:
		mxSetField(mstate,0,"algorithm_class",
		    mxCreateString("LineSearch"));
		break;
	}

	switch(state.Minv_type){
	    case peopt::Identity_t:
		mxSetField(mstate,0,"Minv_type",
		    mxCreateString("Identity"));
		break;
	    case peopt::InvBFGS_t:
		mxSetField(mstate,0,"Minv_type",
		    mxCreateString("InvBFGS"));
		break;
	    case peopt::InvSR1_t:
		mxSetField(mstate,0,"Minv_type",
		    mxCreateString("InvSR1"));
		break;
	    case peopt::External_t:
		mxSetField(mstate,0,"Minv_type",
		    mxCreateString("External"));
		break;
	    default:
	    	pelab::VS::error("Invalid preconditioner type found in the "
		    "defaults.");
	}

	switch(state.H_type){
	    case peopt::Identity_t:
		mxSetField(mstate,0,"H_type",
		    mxCreateString("Identity"));
		break;
	    case peopt::ScaledIdentity_t:
		mxSetField(mstate,0,"H_type",
		    mxCreateString("ScaledIdentity"));
		break;
	    case peopt::BFGS_t:
		mxSetField(mstate,0,"H_type",
		    mxCreateString("BFGS"));
		break;
	    case peopt::SR1_t:
		mxSetField(mstate,0,"H_type",
		    mxCreateString("SR1"));
		break;
	    case peopt::External_t:
		mxSetField(mstate,0,"H_type",
		    mxCreateString("External"));
		break;
	    default:
	    	pelab::VS::error("Invalid preconditioner type found in the "
		    "defaults.");
	}

	mxSetField(mstate,0,"verbose",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"verbose")))[0]
	    =state.verbose;

	mxSetField(mstate,0,"delta",mxCreateDoubleScalar(state.delta)); 
	mxSetField(mstate,0,"delta_max",
	    mxCreateDoubleScalar(state.delta_max)); 
	mxSetField(mstate,0,"eta1",mxCreateDoubleScalar(state.eta1)); 
	mxSetField(mstate,0,"eta2",mxCreateDoubleScalar(state.eta2)); 
	mxSetField(mstate,0,"alpha",mxCreateDoubleScalar(state.alpha)); 
	
	mxSetField(mstate,0,"linesearch_iter_max",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"linesearch_iter_max")))[0]
	    =state.linesearch_iter_max;

	mxSetField(mstate,0,"eps_ls",mxCreateDoubleScalar(state.eps_ls)); 

	switch(state.dir){
	    case peopt::SteepestDescent_t:
		mxSetField(mstate,0,"dir",
		    mxCreateString("SteepestDescent"));
		break;
	    case peopt::FletcherReeves_t:
		mxSetField(mstate,0,"dir",mxCreateString("FletcherReeves"));
		break;
	    case peopt::PolakRibiere_t:
		mxSetField(mstate,0,"dir",mxCreateString("PolakRibiere"));
		break;
	    case peopt::HestenesStiefel_t:
		mxSetField(mstate,0,"dir",
		    mxCreateString("HestenesStiefel"));
		break;
	    case peopt::LimitedMemoryBFGS_t:
		mxSetField(mstate,0,"dir",
		    mxCreateString("LimitedMemoryBFGS"));
		break;
	    case peopt::NewtonCG_t:
		mxSetField(mstate,0,"dir",mxCreateString("NewtonCG"));
		break;
	}

	switch(state.kind){
	    case peopt::Brents_t:
		mxSetField(mstate,0,"kind",mxCreateString("Brents"));
		break;
	    case peopt::GoldenSection_t:
		mxSetField(mstate,0,"kind",mxCreateString("GoldenSection"));
		break;
	    case peopt::BackTracking_t:
		mxSetField(mstate,0,"kind",mxCreateString("BackTracking"));
		break;
	    case peopt::TwoPointA_t:
		mxSetField(mstate,0,"kind",mxCreateString("TwoPointA"));
		break;
	    case peopt::TwoPointB_t:
		mxSetField(mstate,0,"kind",mxCreateString("TwoPointB"));
		break;
	}

	// In the case that we're just doing a problem setup, don't return
	// all of the parameters.  
	if(brief) return; 

	mxSetField(mstate,0,"iter",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"iter")))[0]=state.iter;

	switch(state.opt_stop){
	    case peopt::NotConverged:
		mxSetField(mstate,0,"opt_stop",
		    mxCreateString("NotConverged"));
		break;
	    case peopt::RelativeGradientSmall:
		mxSetField(mstate,0,"opt_stop",
		    mxCreateString("RelativeGradientSmall"));
		break;
	    case peopt::RelativeStepSmall:
		mxSetField(mstate,0,"opt_stop",
		    mxCreateString("RelativeStepSmall"));
		break;
	    case peopt::MaxItersExceeded:
		mxSetField(mstate,0,"opt_stop",
		    mxCreateString("MaxItersExceeded"));
		break;
	}

	mxSetField(mstate,0,"krylov_iter",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"krylov_iter")))[0]
	    =state.krylov_iter;

	mxSetField(mstate,0,"krylov_iter_total",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"krylov_iter_total")))[0]
	    =state.krylov_iter_total;
	
	switch(state.krylov_stop){
	    case peopt::NotConverged:
		mxSetField(mstate,0,"krylov_stop",
		    mxCreateString("NegativeCurvature"));
		break;
	    case peopt::RelativeErrorSmall:
		mxSetField(mstate,0,"krylov_stop",
		    mxCreateString("RelativeErrorSmall"));
		break;
	    case peopt::MaxKrylovItersExceeded:
		mxSetField(mstate,0,"krylov_stop",
		    mxCreateString("MaxKrylovItersExceeded"));
		break;
	    case peopt::TrustRegionViolated:
		mxSetField(mstate,0,"krylov_stop",
		    mxCreateString("TrustRegionViolated"));
		break;
	}

	mxSetField(mstate,0,"krylov_rel_err",
	    mxCreateDoubleScalar(state.krylov_rel_err)); 
	mxSetField(mstate,0,"norm_g",mxCreateDoubleScalar(state.norm_g)); 
	mxSetField(mstate,0,"norm_gtyp",
	    mxCreateDoubleScalar(state.norm_gtyp)); 
	mxSetField(mstate,0,"norm_s",mxCreateDoubleScalar(state.norm_s)); 
	mxSetField(mstate,0,"norm_styp",
	    mxCreateDoubleScalar(state.norm_styp)); 
	mxSetField(mstate,0,"u",mxDuplicateArray(*(state.u.begin())));
	mxSetField(mstate,0,"g",mxDuplicateArray(*(state.g.begin())));
	mxSetField(mstate,0,"s",mxDuplicateArray(*(state.s.begin())));
	mxSetField(mstate,0,"u_old",mxDuplicateArray(*(state.u_old.begin())));
	mxSetField(mstate,0,"g_old",mxDuplicateArray(*(state.g_old.begin())));
	mxSetField(mstate,0,"s_old",mxDuplicateArray(*(state.s_old.begin())));

	mxSetField(mstate,0,"oldY",mxCreateCellMatrix(1,state.oldY.size()));
	VectorList::iterator y=state.oldY.begin();
	for(int i=0;i<state.oldY.size();i++){
	    mxSetCell(mxGetField(mstate,0,"oldY"),i,mxDuplicateArray(*y));
	    y++;
	}
	
	mxSetField(mstate,0,"oldS",mxCreateCellMatrix(1,state.oldS.size()));
	VectorList::iterator s=state.oldS.begin();
	for(int i=0;i<state.oldS.size();i++){
	    mxSetCell(mxGetField(mstate,0,"oldS"),i,mxDuplicateArray(*s));
	    s++;
	}
	
	mxSetField(mstate,0,"obj_u",mxCreateDoubleScalar(state.obj_u)); 
	mxSetField(mstate,0,"obj_ups",mxCreateDoubleScalar(state.obj_ups)); 
	mxSetField(mstate,0,"rho",mxCreateDoubleScalar(state.rho)); 

	mxSetField(mstate,0,"rejected_trustregion",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"rejected_trustregion")))[0]
	    =state.rejected_trustregion;

	mxSetField(mstate,0,"linesearch_iter",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"linesearch_iter")))[0]
	    =state.linesearch_iter;
	
	mxSetField(mstate,0,"linesearch_iter_total",
	    mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL));
	((int*)mxGetPr(mxGetField(mstate,0,"linesearch_iter_total")))[0]
	    =state.linesearch_iter_total;
    }
    
    // Converts a matlab struct to a peopt state
    static void structToState(
	const mxArray* mstate,
	peopt::core <VS>::State& state,
	const bool brief=false
    ){ 

	// Override elements in the state with those provided by the user
	for(int i=0;i<mxGetNumberOfFields(mstate);i++) {

	    // Get the field name
	    const std::string fname=mxGetFieldNameByNumber(mstate,i);

	    // Figure out what property this corresponds to
	    if(fname=="eps_g")
	    	state.eps_g=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="eps_s")
	    	state.eps_s=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="stored_history")
	    	state.stored_history=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="history_reset")
	    	state.history_reset=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="iter_max")
	    	state.iter_max=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="krylov_iter_max")
	    	state.krylov_iter_max=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="eps_krylov")
	    	state.eps_krylov=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="algorithm_class"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="TrustRegion")
			state.algorithm_class=peopt::TrustRegion;
		    else if(value=="LineSearch")
			state.algorithm_class=peopt::LineSearch;
		    else
			pelab::VS::error(
			    "Invalid algorithm class: "+value+"\n");
		}
	    } else if(fname=="Minv_type"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="Identity")
			state.Minv_type=peopt::Identity_t;
		    else if(value=="InvBFGS")
			state.Minv_type=peopt::InvBFGS_t;
		    else if(value=="InvSR1")
			state.Minv_type=peopt::InvSR1_t;
		    else if(value=="External")
			state.Minv_type=peopt::External_t;
		    else
			pelab::VS::error(
			    "Invalid preconditioner: "+value+"\n");
		}
	    } else if(fname=="H_type"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="Identity")
			state.H_type=peopt::Identity_t;
		    else if(value=="ScaledIdentity")
			state.H_type=peopt::ScaledIdentity_t;
		    else if(value=="BFGS")
			state.H_type=peopt::BFGS_t;
		    else if(value=="SR1")
			state.H_type=peopt::SR1_t;
		    else if(value=="External")
			state.H_type=peopt::External_t;
		    else
			pelab::VS::error(
			    "Invalid Hessian: "+value+"\n");
		}
	    }
	    else if(fname=="verbose")
	    	state.verbose=(int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="delta")
	    	state.delta=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="delta_max")
	    	state.delta_max=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="eta1")
	    	state.eta1=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="eta2")
	    	state.eta2=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="alpha")
	    	state.alpha=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="linesearch_iter_max")
	    	state.linesearch_iter_max=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="eps_ls")
	    	state.eps_ls=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="dir"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="SteepestDescent")
			state.dir=peopt::SteepestDescent_t;
		    else if(value=="FletcherReeves")
			state.dir=peopt::FletcherReeves_t;
		    else if(value=="PolakRibiere")
			state.dir=peopt::PolakRibiere_t;
		    else if(value=="HestenesStiefel")
			state.dir=peopt::HestenesStiefel_t;
		    else if(value=="LimitedMemoryBFGS")
			state.dir=peopt::LimitedMemoryBFGS_t;
		    else if(value=="NewtonCG")
			state.dir=peopt::NewtonCG_t;
		    else
			pelab::VS::error(
			    "Invalid linesearch direction: "+value+"\n");
		}
	    } else if(fname=="kind"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="Brents")
			state.kind=peopt::Brents_t;
		    else if(value=="GoldenSection")
			state.kind=peopt::GoldenSection_t;
		    else if(value=="BackTracking")
			state.kind=peopt::BackTracking_t;
		    else if(value=="TwoPointA")
			state.kind=peopt::TwoPointA_t;
		    else if(value=="TwoPointB")
			state.kind=peopt::TwoPointB_t;
		    else
			pelab::VS::error(
			    "Invalid linesearch kind: "+value+"\n");
		}
	    } else if(fname!="F" && fname!="G" && fname !="H" &&
		fname!="Minv" && fname!="StateManipulator" && brief!=false)
	    	pelab::VS::error("Invalid parameter: "+fname+"\n");

	    // If we're just doing initialization, skip the rest of the
	    // parameters
	    if(brief || fname=="eps_g" || fname=="eps_s" ||
		fname=="stored_history" || fname=="history_reset" ||
		fname=="iter_max" || fname=="krylov_iter_max" ||
		fname=="eps_krylov" || fname=="algorithm_class" ||
		fname=="Minv_type" || fname=="H_type" || fname=="verbose" ||
		fname=="delta" || fname=="delta_max" ||
		fname=="eta1" || fname=="eta2" ||
		fname=="alpha" || fname=="linesearch_iter_max" ||
		fname=="eps_ls" || fname=="dir" ||
		fname=="kind") continue;
	    
	    if(fname=="iter")
	    	state.iter=(int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="opt_stop"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="NotConverged")
			state.opt_stop=peopt::NotConverged;
		    else if(value=="RelativeGradientSmall")
			state.opt_stop=peopt::RelativeGradientSmall;
		    else if(value=="RelativeStepSmall")
			state.opt_stop=peopt::RelativeStepSmall;
		    else if(value=="MaxItersExceeded")
			state.opt_stop=peopt::MaxItersExceeded;
		    else
			pelab::VS::error(
			    "Invalid optimization stopping condition: "+
			    value+"\n");
		}
	    } else if(fname=="krylov_iter")
	    	state.krylov_iter=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="krylov_iter_total")
	    	state.krylov_iter_total=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="krylov_stop"){
	    	mxArray* field=mxGetFieldByNumber(mstate,0,i);
		if(field){
		    int buflen = mxGetN(field);
		    std::string value(buflen,'\0');
		    mxGetString(field,&(value[0]),buflen+1);
		    if(value=="NegativeCurvature")
			state.krylov_stop=peopt::NegativeCurvature;
		    else if(value=="RelativeErrorSmall")
			state.krylov_stop=peopt::RelativeErrorSmall;
		    else if(value=="MaxKrylovItersExceeded")
			state.krylov_stop=peopt::MaxKrylovItersExceeded;
		    else if(value=="TrustRegionViolated")
			state.krylov_stop=peopt::TrustRegionViolated;
		    else
			pelab::VS::error(
			    "Invalid Krylov method stopping condition: "+
			    value+"\n");
		}
	    } else if(fname=="krylov_rel_err")
	    	state.krylov_rel_err=
		    mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="norm_g")
	    	state.norm_g=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="norm_gtyp")
	    	state.norm_gtyp=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="norm_s")
	    	state.norm_s=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="norm_styp")
	    	state.norm_styp=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="u")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.u.begin()));
	    else if(fname=="g")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.g.begin()));
	    else if(fname=="s")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.s.begin()));
	    else if(fname=="u_old")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.u_old.begin()));
	    else if(fname=="g_old")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.g_old.begin()));
	    else if(fname=="s_old")
	    	VS::copy(mxGetFieldByNumber(mstate,0,i),*(state.s_old.begin()));
	    else if(fname=="oldY") {
		VectorList::iterator y=state.oldY.begin();
		for(int i=0;i<state.oldY.size();i++){
		    VS::copy(mxGetCell(mxGetField(mstate,0,"oldY"),i),*y);
		    y++;
		}
	    } else if(fname=="oldS") {
		VectorList::iterator s=state.oldS.begin();
		for(int i=0;i<state.oldS.size();i++){
		    VS::copy(mxGetCell(mxGetField(mstate,0,"oldS"),i),*s);
		    s++;
		}
	    } else if(fname=="obj_u")
	    	state.obj_u=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="obj_ups")
	    	state.obj_ups=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="rho")
	    	state.rho=mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="rejected_trustregion")
	    	state.rejected_trustregion=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="linesearch_iter")
	    	state.linesearch_iter=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else if(fname=="linesearch_iter_total")
	    	state.linesearch_iter_total=
		    (int)mxGetScalar(mxGetFieldByNumber(mstate,0,i));
	    else 
	    	pelab::VS::error("Invalid parameter: "+fname+"\n");
	}
    }

    // Implements a differentiable operator using Matlab functions
    class DiffOperator: public peopt::DiffOperator <VS,VS> {
    private:
    	mxArray*& g;
	mxArray*& gp;
	mxArray*& gps;
	mxArray*& gpps;
	const int mindex;
    public:
    	DiffOperator(
	    mxArray*& g_,
	    mxArray*& gp_,
	    mxArray*& gps_,
	    mxArray*& gpps_,
	    const int mindex_)
	    : g(g_), gp(gp_), gps(gps_), gpps(gpps_), mindex(mindex_) {}

	// Basic application
	class F : public peopt::Operator <VS,VS> {
	private:
	    mxArray*& f;
	    int i;
	public:
	    F(mxArray*& f_,int i_) : f(f_), i(i_) {};	
	    void operator () (const VS::Vector& x,VS::Vector &y) const {
		
		// Create memory for the index 
		mxArray* idx;
		idx=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
		((int*)mxGetPr(idx))[0]=i+1;

		// Create the inputs 
		mxArray* input[3]={f,idx,x};
		
		// Compute the gradient function
		mexCallMATLAB(1,&y,3,input,"feval");
		
		// Free memory for the index
		mxDestroyArray(idx);
	    }
	};
	std::auto_ptr < peopt::Operator <VS,VS> > f(const int i) const {
	    return std::auto_ptr < peopt::Operator <VS,VS> > (new F(g,i));
	}

	// Derivative of the operator 
	class Fp : public peopt::Operator <VS,VS> {
	private:
	    mxArray*& fp;
	    const VS::Vector& x;
	    int i;
	public:
	    Fp(mxArray*& fp_,int i_,const VS::Vector& x_)
		: fp(fp_), i(i_), x(x_) {};	
	    void operator () (const VS::Vector& eta,VS::Vector &y) const {
		
		// Create memory for the index 
		mxArray* idx;
		idx=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
		((int*)mxGetPr(idx))[0]=i+1;

		// Create the inputs 
		mxArray* input[4]={fp,idx,x,eta};
		
		// Compute the gradient function
		mexCallMATLAB(1,&y,4,input,"feval");
		
		// Free memory for the index
		mxDestroyArray(idx);
	    }
	};
	std::auto_ptr < peopt::Operator <VS,VS> >
	    fp(const int i, const VS::Vector& x
	) const {
	    return std::auto_ptr < peopt::Operator <VS,VS> > (new Fp(gp,i,x));
	}

	// Derivative adjoint 
	class Fps : public peopt::Operator <VS,VS> {
	private:
	    mxArray*& fps;
	    const VS::Vector& x;
	    int i;
	public:
	    Fps(mxArray*& fps_,int i_,const VS::Vector& x_)
		: fps(fps_), i(i_), x(x_) {};	
	    void operator () (const VS::Vector& xi,VS::Vector &y) const {

		// Create memory for the index 
		mxArray* idx;
		idx=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
		((int*)mxGetPr(idx))[0]=i+1;

		// Create the inputs 
		mxArray* input[4]={fps,idx,x,xi};
		
		// Compute the gradient function
		mexCallMATLAB(1,&y,4,input,"feval");
		
		// Free memory for the index
		mxDestroyArray(idx);
	    }
	};
	std::auto_ptr < peopt::Operator <VS,VS> >
	    fps(const int i, const VS::Vector& x
	) const {
	    return std::auto_ptr < peopt::Operator <VS,VS> >(new Fps(gps,i,x));
	}
	
	// Second derivative in the direction eta, adjoint, in the direction xi
	class Fpps : public peopt::Operator <VS,VS> {
	private:
	    mxArray*& fpps;
	    const VS::Vector& x;
	    const VS::Vector& xi;
	    int i;
	public:
	    Fpps(
		mxArray*& fpps_,	
		int i_,
		const VS::Vector& x_,
		const VS::Vector& xi_)
		: fpps(fpps_), i(i_), x(x_), xi(xi_) {};	
	    void operator () (const VS::Vector& eta,VS::Vector &y) const {

		// Create memory for the index 
		mxArray* idx;
		idx=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
		((int*)mxGetPr(idx))[0]=i+1;

		// Create the inputs 
		mxArray* input[5]={fpps,idx,x,eta,xi};
		
		// Compute the gradient function
		mexCallMATLAB(1,&y,5,input,"feval");
		
		// Free memory for the index
		mxDestroyArray(idx);
	    }
	};
	std::auto_ptr < peopt::Operator <VS,VS> >
	    fpps(const int i, const VS::Vector& x, const VS::Vector& xi
	) const {
	    return std::auto_ptr < peopt::Operator <VS,VS> >
		(new Fpps(gpps,i,x,xi));
	}

	// Returns the maximum index.  In other words, the total number of
	// PDEs under consideration.
	int max_index() const { return mindex; }
    };
};

const std::string pelab::VS::copy_="copy";
const std::string pelab::VS::scal_="scal";
const std::string pelab::VS::axpy_="axpy";
const std::string pelab::VS::innr_="innr";
const std::string pelab::VS::init_="init";


#endif
