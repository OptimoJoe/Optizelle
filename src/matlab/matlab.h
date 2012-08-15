#ifndef matlab_h
#define matlab_h

#include <string>
#include <sstream>
#include "mex.h"
#include "peopt/peopt.h"

// Defines the vector space used for optimization.
#define CREATE_VS(name) \
template <typename Real> \
struct name { \
    typedef mxArray* Vector; \
\
    static mxArray* init_; \
    static mxArray* copy_; \
    static mxArray* scal_; \
    static mxArray* zero_; \
    static mxArray* axpy_; \
    static mxArray* innr_; \
    static mxArray* prod_; \
    static mxArray* id_; \
    static mxArray* linv_; \
    static mxArray* barr_; \
    static mxArray* srch_; \
\
    /* Memory allocation and size setting */ \
    static void init(const Vector& x, Vector& y) { \
        Vector input[2]={init_,x}; \
        mexCallMATLAB(1,&y,2,input,"feval"); \
    } \
\
    /* y <- x (Shallow.  No memory allocation.) */ \
    static void copy(const Vector& x, Vector& y) { \
        Vector input[2]={copy_,x}; \
        mexCallMATLAB(1,&y,2,input,"feval"); \
    } \
\
    /* x <- alpha * x */ \
    static void scal(const Real& alpha_, Vector& x) { \
        /* Create memory for the scalar */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
        mxGetPr(alpha)[0]=alpha_; \
\
        /* Create the inputs and outputs */ \
        Vector input[3]={scal_,alpha,x}; \
\
        /* Compute the scalar multiplication */ \
        mexCallMATLAB(1,&x,3,input,"feval"); \
\
        /* Free memory from the scalar */ \
        mxDestroyArray(alpha); \
    } \
\
    /* x <- 0 */ \
    static void zero(Vector& x) { \
        Vector input[2]={zero_,x}; \
        mexCallMATLAB(1,&x,2,input,"feval"); \
    } \
\
    /* y <- alpha * x + y */ \
    static void axpy(const Real& alpha_, const Vector& x, Vector& y) { \
        /* Create memory for the scalar */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
        mxGetPr(alpha)[0]=alpha_; \
\
        /* Create the inputs and outputs */ \
        Vector input[4]={axpy_,alpha,x,y}; \
\
        /* Compute the addition */ \
        mexCallMATLAB(1,&y,4,input,"feval"); \
\
        /* Free memory from the scalar */ \
        mxDestroyArray(alpha); \
    } \
\
    /* innr <- <x,y> */ \
    static Real innr(const Vector& x,const Vector& y) { \
        /* Create memory for the result */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
\
        /* Create the inputs and outputs */ \
        Vector input[3]={innr_,x,y}; \
\
        /* Compute the inner product */ \
        mexCallMATLAB(1,&alpha,3,input,"feval"); \
\
        /* Get the result of the computation */ \
        double alpha_=mxGetScalar(alpha); \
\
        /* Free memory from the scalar */ \
        mxDestroyArray(alpha); \
\
        /* Return the result */ \
        return alpha_; \
    } \
\
    /* Jordan product, z <- x o y */ \
    static void prod(const Vector& x, const Vector& y, Vector& z) { \
        /* Create the inputs and outputs */ \
        Vector input[3]={prod_,x,y}; \
\
        /* Compute the product */ \
        mexCallMATLAB(1,&z,3,input,"feval"); \
    } \
\
    /* Identity element, x <- e such that x o e = x */ \
    static void id(Vector& x) { \
        Vector input[2]={id_,x}; \
        mexCallMATLAB(1,&x,2,input,"feval"); \
    } \
\
    /* Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y */ \
    static void linv(const Vector& x, const Vector& y, Vector& z) { \
        /* Create the inputs and outputs */ \
        Vector input[3]={linv_,x,y}; \
\
        /* Compute the product */ \
        mexCallMATLAB(1,&z,3,input,"feval"); \
    } \
\
    /* Barrier function, barr <- barr(x) where x o grad barr(x) = e */ \
    static Real barr(const Vector& x,const Vector& y) { \
        /* Create memory for the result */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
\
        /* Create the inputs and outputs */ \
        Vector input[3]={barr_,x,y}; \
\
        /* Compute the barrier function */ \
        mexCallMATLAB(1,&alpha,3,input,"feval"); \
\
        /* Get the result of the computation */ \
        double alpha_=mxGetScalar(alpha); \
\
        /* Free memory from the scalar */ \
        mxDestroyArray(alpha); \
\
        /* Return the result */ \
        return alpha_; \
    } \
\
    /* Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0} */ \
    /* where y > 0.  If the argmax is infinity, then return Real(-1.). */ \
    static Real srch(const Vector& x,const Vector& y) { \
        /* Create memory for the result */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
\
        /* Create the inputs and outputs */ \
        Vector input[3]={srch_,x,y}; \
\
        /* Compute the search function */ \
        mexCallMATLAB(1,&alpha,3,input,"feval"); \
\
        /* Get the result of the computation */ \
        double alpha_=mxGetScalar(alpha); \
\
        /* Free memory from the scalar */ \
        mxDestroyArray(alpha); \
\
        /* Return the result */ \
        return alpha_; \
    } \
}; \
template <typename Real> mxArray* name<Real>::init_; \
template <typename Real> mxArray* name<Real>::copy_; \
template <typename Real> mxArray* name<Real>::scal_; \
template <typename Real> mxArray* name<Real>::zero_; \
template <typename Real> mxArray* name<Real>::axpy_; \
template <typename Real> mxArray* name<Real>::innr_; \
template <typename Real> mxArray* name<Real>::prod_; \
template <typename Real> mxArray* name<Real>::id_; \
template <typename Real> mxArray* name<Real>::linv_; \
template <typename Real> mxArray* name<Real>::barr_; \
template <typename Real> mxArray* name<Real>::srch_;

#define SET_VS(name,from) \
    name <double>::init_=mxGetField(from,0,"init"); \
    name <double>::copy_=mxGetField(from,0,"copy"); \
    name <double>::scal_=mxGetField(from,0,"scal"); \
    name <double>::zero_=mxGetField(from,0,"zero"); \
    name <double>::axpy_=mxGetField(from,0,"axpy"); \
    name <double>::innr_=mxGetField(from,0,"innr"); 

#define SET_EJA(name,from) \
    name <double>::init_=mxGetField(from,0,"prod"); \
    name <double>::copy_=mxGetField(from,0,"id"); \
    name <double>::scal_=mxGetField(from,0,"linv"); \
    name <double>::zero_=mxGetField(from,0,"barr"); \
    name <double>::axpy_=mxGetField(from,0,"srch"); 


// A messaging utility that hooks directly into Matlab
struct MatlabMessaging : public peopt::Messaging {
    // Prints a message
    void print(const std::string msg) const {
        std::string msg_ = msg + "\n";
        mexPrintf(msg_.c_str());
    }

    // Prints an error
    void error(const std::string msg) const {
        mexErrMsgTxt(msg.c_str());
    }
};

// Checks that a structure array contains certain functions
void check_fns(
    const mxArray* X,
    std::string ops[],
    int nops,
    std::string type,
    std::string name
) {
    // The argument should be a structure array
    if(!mxIsStruct(X))
        MatlabMessaging().error("The " + type + " " + name +
            " must be a structure array.");

    // The vector space must contain the elements:
    // init, copy, scal, zero, axpy, innr
    for(int i=0;i<nops;i++)
        if(mxGetField(X,0,ops[i].c_str())==NULL)
            MatlabMessaging().error("Missing the " + ops[i] + 
                " function in the " + type + " " + name +".");
    
    // Check that all of the vector space functions are really functions
    for(int i=0;i<nops;i++)
        if(!mxIsClass(mxGetField(X,0,ops[i].c_str()),"function_handle"))
            MatlabMessaging().error("The field " + ops[i] +
                " in the " + type + " " + name + " must be a function.");
}

// Check that we have all the operations necessary for a vector space
void check_vector_space(const mxArray* X,std::string name) {
    std::string ops[6] = { "init", "copy", "scal", "zero", "axpy", "innr" };
    check_fns(X,ops,6,"vector space",name);
}

// Check that we have all the operations necessary for a Euclidean-Jordan
// algebra
void check_euclidean_jordan(const mxArray* X,std::string name) {
    // First, check that we have a valid vector space
    check_vector_space(X,name);

    // Next, we check we have a valid Euclidean-Jordan algebra
    std::string ops[5] = { "prod", "id", "linv", "barr", "srch" };
    check_fns(X,ops,5,"vector space",name);
}

// Check that we have all the operations necessary for a scalar valued function 
void check_scalar_valued_fn(const mxArray* X,std::string name) {
    std::string ops[3] = { "eval", "grad", "hessvec" };
    check_fns(X,ops,3,"scalar valued function",name);
}


// Check that we have all the operations necessary for a vector valued function 
void check_vector_valued_fn(const mxArray* X,std::string name) {
    std::string ops[4] = { "eval", "eval_p", "eval_ps", "eval_pps" };
    check_fns(X,ops,4,"vector valued function",name);
}

// Determine the problem class from a bundle of vector spaces 
peopt::ProblemClass::t get_problem_class_from_vs(const mxArray* VS) {
    // First, check that we have a structure array        
    if(!mxIsStruct(VS))
        MatlabMessaging().error("The bundle of vector spaces must be a "
            "structure array.");

    // Next, insure that we have at least X
    if(mxGetField(VS,0,"X")==NULL)
        MatlabMessaging().error("The bundle of vector spaces must " 
            "contain X.");

    // Check that X is a vector space
    check_vector_space(mxGetField(VS,0,"X"),"X");

    // Set the initial algorithm class
    peopt::ProblemClass::t problem_class
        = peopt::ProblemClass::Unconstrained;

    // Next, check if we have any equality constraints
    if(mxGetField(VS,0,"Y") != NULL) {

        // Check for the vector space operations
        check_vector_space(mxGetField(VS,0,"Y"),"Y");

        // Adjust the algorithm class
        problem_class
            = peopt::ProblemClass::EqualityConstrained;
    }

    // Next, check if we have inequality constraints
    if(mxGetField(VS,0,"Z") != NULL) {

        // Check for the Euclidean-Jordan algebra operations
        check_euclidean_jordan(mxGetField(VS,0,"Z"),"Z");

        // Adjust the algorithm class
        if(problem_class==peopt::ProblemClass::EqualityConstrained)
            problem_class = peopt::ProblemClass::Constrained;
        else
            problem_class = peopt::ProblemClass::InequalityConstrained;
    }

    return problem_class;
}

// Determine the problem class from a bundle of vector spaces 
peopt::ProblemClass::t get_problem_class_from_fns(const mxArray* fns) {
    // First, check that we have a structure array        
    if(!mxIsStruct(fns))
        MatlabMessaging().error("The bundle of functions must be a "
            "structure array.");

    // Next, insure that we have at least f
    if(mxGetField(fns,0,"f")==NULL)
        MatlabMessaging().error("The bundle of functions must " 
            "contain f.");

    // Check that f is a scalar value function 
    check_scalar_valued_fn(mxGetField(fns,0,"f"),"f");

    // Set the initial algorithm class
    peopt::ProblemClass::t problem_class
        = peopt::ProblemClass::Unconstrained;

    // Next, check if we have any equality constraints
    if(mxGetField(fns,0,"g") != NULL) {

        // Check for the vector space operations
        check_vector_valued_fn(mxGetField(fns,0,"g"),"g");

        // Adjust the algorithm class
        problem_class
            = peopt::ProblemClass::EqualityConstrained;
    }

    // Next, check if we have inequality constraints
    if(mxGetField(fns,0,"h") != NULL) {

        // Check for the Euclidean-Jordan algebra operations
        check_vector_valued_fn(mxGetField(fns,0,"h"),"h");

        // Adjust the algorithm class
        if(problem_class==peopt::ProblemClass::EqualityConstrained)
            problem_class = peopt::ProblemClass::Constrained;
        else
            problem_class = peopt::ProblemClass::InequalityConstrained;
    }

    return problem_class;
}


#endif
