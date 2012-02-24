#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage="\nParameter estimation Newton Hessian-vector"
	"product: getNewton(r,h,vecs,u,du)\n";

    // Check inputs and outputs 
    if(nInput!=5)
        mexErrMsgTxt((usage+"Incorrect number of input arguments.\n").c_str());
    if(nOutput!=1)
        mexErrMsgTxt((usage+"Incorrect number of output arguments.\n").c_str());

    // Get the residual and solution operators
    mxArray* R=const_cast<mxArray*>(pInput[0]);
    	mxArray* r_=mxGetField(R,0,"f");
    	mxArray* rp=mxGetField(R,0,"fp");
    	mxArray* rps=mxGetField(R,0,"fps");
    	mxArray* rpps=mxGetField(R,0,"fpps");
	int rmax=int((mxGetPr(mxGetField(R,0,"max_index")))[0]);
    mxArray* H=const_cast<mxArray*>(pInput[1]);
    	mxArray* h_=mxGetField(H,0,"f");
    	mxArray* hp=mxGetField(H,0,"fp");
    	mxArray* hps=mxGetField(H,0,"fps");
    	mxArray* hpps=mxGetField(H,0,"fpps");
	int hmax=int((mxGetPr(mxGetField(H,0,"max_index")))[0]);

    // Get a representative sample of vectors.  Then, extract each vector.
    mxArray* vecs=const_cast<mxArray*>(pInput[2]);
    mxArray* vec_y=mxGetField(vecs,0,"y");
    mxArray* vec_z=mxGetField(vecs,0,"z");

    // Get the point where we're evaluating the objective
    mxArray* u=const_cast<mxArray*>(pInput[3]);

    // Get the direction in which we're taking the Hessian-vector product
    mxArray* eta=const_cast<mxArray*>(pInput[4]);

    // Create the residual and solution operators
    pelab::DiffOperator r(r_,rp,rps,rpps,rmax);
    pelab::DiffOperator h(h_,hp,hps,hpps,hmax);

    // Now, create the objective function based on these operators
    peopt::parest<pelab::VS,pelab::VS,pelab::VS>::Newton
	HH(r,h,u,vec_y,vec_z);

    // Create memory for the hess-vec product 
    mxArray* hv=mxDuplicateArray(u);

    // Evaluate the hess-vec 
    HH(eta,hv);

    // Return the result to matlab
    pOutput[0]=hv;
}
