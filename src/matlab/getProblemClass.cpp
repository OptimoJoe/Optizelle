#include <string>
#include <sstream>
#include "mex.h"
#include "matlab.h"

CREATE_VS(XX);
CREATE_VS(YY);

template <typename Real, template <typename> class XX>
void someFn(
    const typename XX <Real>::Vector& x,
    typename XX<Real>::Vector& y
) {
    typedef XX <Real> X;
    X::copy(x,y);
    X::scal(2.,y);
    X::axpy(1.,x,y);
    std::stringstream msg;
    msg << "This inner product is: " <<  X::innr(y,y) << std::endl;
    mexPrintf(msg.str().c_str());
}

void mexFunction(
    int nOutput,mxArray* pOutput[],
    int nInput,const mxArray* pInput[]
) {
    // Argument Checking:

    // Check the number of arguments
    if(nInput!=2 || nOutput!=0)
        MatlabMessaging().error("Invalid number of arguments.");

    // Get the problem class
    peopt::ProblemClass::t problem_class_vs
        =get_problem_class_from_vs(pInput[0]);
    peopt::ProblemClass::t problem_class_fn
        =get_problem_class_from_fns(pInput[1]);

    // Setup the vector space
    //SET_VS(XX,X);

    // Print out the problem class
    MatlabMessaging().print(peopt::ProblemClass::to_string(problem_class_vs));
    MatlabMessaging().print(peopt::ProblemClass::to_string(problem_class_fn));
}
