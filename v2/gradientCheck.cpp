#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage="\nFinite difference gradient check: "
	"gradientCheck(f,g,u,eta)\n";

    // Check inputs and outputs 
    if(nInput!=4)
        mexErrMsgTxt((usage+"Incorrect number of input arguments.\n").c_str());
    if(nOutput!=0)
        mexErrMsgTxt((usage+"Incorrect number of output arguments.\n").c_str());

    // Create a shortcut for the base of the finite difference check 
    mxArray* u=const_cast<mxArray*>(pInput[2]);

    // Create a shortcut for the direction
    mxArray* eta=const_cast<mxArray*>(pInput[3]);

    // Next, get the objective and gradient
    mxArray* f_=const_cast <mxArray*>(pInput[0]);
    mxArray* g_=const_cast <mxArray*>(pInput[1]);
    
    // Now, create a functional and gradient based on this information
    pelab::Functional F(f_);
    pelab::Gradient G(g_);

    // Perform the finite difference check 
    pelab::VS::print("Finite difference check on the gradient of the "
	"objective.");
    peopt::derivativeCheck<pelab::VS>(F,G,u,eta);
}
