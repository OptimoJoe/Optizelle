#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage=
	"\nParameter estimation laboratory\n\n"
	"Determine default parameters: params=pelab()\n"
	"Optimize: [y why_stop]=pelab(params,u): y <- min f(u)\n\n";

    // Check inputs and outputs 
    if(nInput!=0 && nInput!=2)
        mexErrMsgTxt((usage+"Incorrect number of input arguments.\n").c_str());
    if(nOutput!=1 && nOutput!=2)
        mexErrMsgTxt((usage+"Incorrect number of output arguments.\n").c_str());

    // In the case that we receive no input parameters, return the current
    // set of default parameters
    if(nInput==0) {
    	peopt::core<pelab::VS>::State state;
	pelab::stateToStruct(state,pOutput[0],true);
    } else {

	// Create a shortcut for the initial guess
	mxArray* u=const_cast<mxArray**>(pInput)[1];
	
	// Set the optimization state
	peopt::core<pelab::VS>::State state(u);

	// Initialize memory for the name of the function, gradient, Hessian,
	// and preconditioner.
	mxArray* F_=NULL;
	mxArray* G_=NULL;
	mxArray* H_=NULL;
	mxArray* Minv_=NULL;

	// Also, keep track of the state manipulator
	mxArray* smanip_=NULL;

	// Setup the initial state
	pelab::structToState(pInput[0],state,true);

	// Find the function, gradient, Hessian, and state manipulator 
	for(int i=0;i<mxGetNumberOfFields(pInput[0]);i++) {

	    // Get the field name
	    const std::string fname=mxGetFieldNameByNumber(pInput[0],i);

	    if(fname=="F"){
	    	F_=mxGetFieldByNumber(pInput[0],0,i);
	    } else if(fname=="G"){
	    	G_=mxGetFieldByNumber(pInput[0],0,i);
	    } else if(fname=="H"){
	    	H_=mxGetFieldByNumber(pInput[0],0,i);
	    } else if(fname=="Minv"){
	    	Minv_=mxGetFieldByNumber(pInput[0],0,i);
	    } else if(fname=="StateManipulator"){
	    	smanip_=mxGetFieldByNumber(pInput[0],0,i);
	    } 
	}

	// Next, insure that we have a function and gradient
	if(!F_)
	    pelab::VS::error(usage+"Objective function not defined.\n");
	if(!G_)
	    pelab::VS::error(usage+"Gradient not defined.\n");

	// Now, we define the objective function and gradient.
	pelab::Functional F(F_);
	pelab::Gradient G(G_);

	// Next, we minimize
	if(state.H_type!=peopt::External_t) 
	    if(!smanip_)
		peopt::core<pelab::VS>::getMin(state,F,G);
	    else {
		pelab::StateManipulator smanip(smanip_);
		peopt::core<pelab::VS>::getMin(state,smanip,F,G);
	    }
	else if(state.Minv_type!=peopt::External_t){
	    pelab::Operator H(H_,*(state.u.begin()));
	    if(!smanip_)
		peopt::core<pelab::VS>::getMin(state,F,G,H);
	    else {
		pelab::StateManipulator smanip(smanip_);
		peopt::core<pelab::VS>::getMin(state,smanip,F,G,H);
	    }
	} else {
	    pelab::Operator H(H_,*(state.u.begin()));
	    pelab::Operator Minv(Minv_,*(state.u.begin()));
	    if(!smanip_) 
		peopt::core<pelab::VS>::getMin(state,F,G,H,Minv);
	    else {
		pelab::StateManipulator smanip(smanip_);
		peopt::core<pelab::VS>::getMin(state,smanip,F,G,H,Minv);
	    }
	}
	
	// Now, copy the result of optimization into the output 
	pOutput[0]=mxDuplicateArray(*(state.u.begin()));

	// If the user requests it, return why the optimization terminated
	if(nOutput==2){
	    switch(state.opt_stop){
	    case peopt::NotConverged:
		pOutput[1]=mxCreateString("NotConverged");
		break;
	    case peopt::RelativeGradientSmall:
		pOutput[1]=mxCreateString("RelativeGradientSmall");
		break;
	    case peopt::RelativeStepSmall:
		pOutput[1]=mxCreateString("RelativeStepSmall");
		break;
	    case peopt::MaxItersExceeded:
		pOutput[1]=mxCreateString("MaxItersExceeded");
		break;
	    }
	}
    }
}
