#include "optizelle.h"

void mexFunction(
    int nOutput,mxArray* pOutput[],
    int nInput,mxArray const * pInput[]
) {
    Optizelle::Matlab::EqualityConstrained::State::readJson(
        nOutput,pOutput,nInput,const_cast <mxArray **> (pInput));
}
