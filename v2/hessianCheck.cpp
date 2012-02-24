#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage="\nFinite difference Hessian check: "
	"hessianCheck(g,H,u,eta)\n";

    // Check inputs and outputs 
    if(nInput!=4)
        mexErrMsgTxt((usage+"Incorrect number of input arguments.\n").c_str());
    if(nOutput!=0)
        mexErrMsgTxt((usage+"Incorrect number of output arguments.\n").c_str());

    // Create a shortcut for the base of the finite difference check 
    mxArray* u=const_cast<mxArray*>(pInput[2]);

    // Create a shortcut for the direction
    mxArray* eta=const_cast<mxArray*>(pInput[3]);

    // Next, get the function names for the gradient and Hessian
    mxArray* g_=const_cast <mxArray*>(pInput[0]);
    mxArray* H_=const_cast <mxArray*>(pInput[1]);
    
    // Now, create a gradient and Hessian based on this information
    pelab::Gradient G(g_);
    pelab::Operator H(H_,u);

    // Perform the finite difference check 
    pelab::VS::print("Finite difference check of the Hessian-vector "
	"product.\n");
    peopt::derivativeCheck<pelab::VS,pelab::VS>(G,H,u,eta,u);
}
