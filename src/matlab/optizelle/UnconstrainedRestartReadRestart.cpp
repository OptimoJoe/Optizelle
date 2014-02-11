#include "optizelle.h"

void mexFunction(
    int nOutput,mxArray* pOutput[],
    int nInput,mxArray const * pInput[]
) {
    Optizelle::Matlab::Unconstrained::Restart::read_restart(
        nOutput,pOutput,nInput,const_cast <mxArray **> (pInput));
}
