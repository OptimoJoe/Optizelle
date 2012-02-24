#include "pelab.h"
#include "mex.h"
#include <string>
#include <iostream>
    
void mexFunction(
    int nOutput, mxArray* pOutput[],
    int nInput, const mxArray* pInput[]
){

    const std::string usage="\nParameter estimation objective value: "
	"getObjective(r,h,vecs,u)\n";

    // Check inputs and outputs 
    if(nInput!=4)
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

    // Create the residual and solution operators
    pelab::DiffOperator r(r_,rp,rps,rpps,rmax);
    pelab::DiffOperator h(h_,hp,hps,hpps,hmax);

    // Now, create the objective function based on these operators
    peopt::parest<pelab::VS,pelab::VS,pelab::VS>::getObjective 
	F(r,h,vec_y,vec_z);

    // Evaluate the objective at the given point
    double alpha=F(u);

    // Return the result to matlab
    pOutput[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetPr(pOutput[0])[0]=alpha;
}
