#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
#include <sstream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage="\nParameter estimation differential operator "
	"check: fpCheck(op,x,eta,xi)\n";

    // Check inputs and outputs 
    if(nInput!=4)
        mexErrMsgTxt((usage+"Incorrect number of input arguments.\n").c_str());
    if(nOutput!=0)
        mexErrMsgTxt((usage+"Incorrect number of output arguments.\n").c_str());

    // Get the residual and solution operators
    mxArray* op=const_cast<mxArray*>(pInput[0]);
    	mxArray* f=mxGetField(op,0,"f");
    	mxArray* fp=mxGetField(op,0,"fp");
    	mxArray* fps=mxGetField(op,0,"fps");
    	mxArray* fpps=mxGetField(op,0,"fpps");
	int fmax=int((mxGetPr(mxGetField(op,0,"max_index")))[0]);

    // Get the point where we're evaluating the objective
    mxArray* x=const_cast<mxArray*>(pInput[1]);

    // Get the direction for the finite difference test
    mxArray* eta=const_cast<mxArray*>(pInput[2]);
    
    // Get some element in the codomain of the operator 
    mxArray* xi=const_cast<mxArray*>(pInput[3]);

    // Create the differentiable operator 
    pelab::DiffOperator Op(f,fp,fps,fpps,fmax);

    // Form the objective function and the gradient.  Then, do a finite
    // difference check on the gradient.  This checks fps. 

    // Do a finite difference check on fp.
    pelab::VS::print("Finite difference check of fp.\n");
    for(int i=0;i<Op.max_index();i++){
    	std::stringstream ss;
	ss << "Calculating " << i+1 << "/" << Op.max_index() << std::endl;
	pelab::VS::print(ss.str());
	peopt::derivativeCheck <pelab::VS,pelab::VS>
	    (*(Op.f(i)),*(Op.fp(i,x)),x,eta,xi);
    }
}
