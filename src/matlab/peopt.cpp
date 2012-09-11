#include <string>
#include <sstream>
#include "mex.h"
#include "peopt/peopt.h"
#include "peopt/json.h"

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
        y=mxDuplicateArray(x); \
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
    static Real barr(const Vector& x) { \
        /* Create memory for the result */ \
        mxArray* alpha; \
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); \
\
        /* Create the inputs and outputs */ \
        Vector input[2]={barr_,x}; \
\
        /* Compute the barrier function */ \
        mexCallMATLAB(1,&alpha,2,input,"feval"); \
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
    name <double>::prod_=mxGetField(from,0,"prod"); \
    name <double>::id_=mxGetField(from,0,"id"); \
    name <double>::linv_=mxGetField(from,0,"linv"); \
    name <double>::barr_=mxGetField(from,0,"barr"); \
    name <double>::srch_=mxGetField(from,0,"srch"); 

CREATE_VS(XX);
CREATE_VS(YY);
CREATE_VS(ZZ);

// A simple scalar valued function using Matlab provided functions 
template <typename Real,template <typename> class XX>
struct MatlabScalarValuedFunction : peopt::ScalarValuedFunction <Real,XX> {
private:
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector Vector;
   
    // The Matlab functions
    mxArray* eval_;
    mxArray* grad_;
    mxArray* hessvec_;

public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabScalarValuedFunction(mxArray* f) {
        eval_=mxGetField(f,0,"eval"); \
        grad_=mxGetField(f,0,"grad"); \
        hessvec_=mxGetField(f,0,"hessvec"); \
    }

    // <- f(x) 
    Real operator () (const Vector& x) const {
        /* Create memory for the result */ 
        mxArray* alpha; 
        alpha=mxCreateDoubleMatrix(1,1,mxREAL); 

        /* Create the inputs and outputs */ 
        Vector input[2]={eval_,x}; 

        /* Compute the inner product */ 
        mexCallMATLAB(1,&alpha,2,input,"feval"); 

        /* Get the result of the computation */ 
        double alpha_=mxGetScalar(alpha); 

        /* Free memory from the scalar */ 
        mxDestroyArray(alpha); 

        /* Return the result */ 
        return alpha_; 
    }

    // g = grad f(x) 
    void grad(const Vector& x,Vector& g) const { 
        Vector input[2]={grad_,x}; 
        mexCallMATLAB(1,&g,2,input,"feval");
    }

    // H_dx = hess f(x) dx 
    void hessvec(const Vector& x,const Vector& dx,Vector& H_dx) const {
        Vector input[3]={hessvec_,x,dx}; 
        mexCallMATLAB(1,&H_dx,3,input,"feval");
    }
};

// A vector valued function using Matlab provided functions 
template <
    typename Real,
    template <typename> class XX,
    template <typename> class YY 
>
struct MatlabVectorValuedFunction : peopt::VectorValuedFunction <Real,XX,YY> {
private:
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector; 
    typedef YY <Real> Y;
    typedef typename Y::Vector Y_Vector; 
   
    // The Matlab functions
    mxArray* eval_;
    mxArray* eval_p_;
    mxArray* eval_ps_;
    mxArray* eval_pps_;

public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabVectorValuedFunction(mxArray* f) {
        eval_=mxGetField(f,0,"eval"); \
        eval_p_=mxGetField(f,0,"eval_p"); \
        eval_ps_=mxGetField(f,0,"eval_ps"); \
        eval_pps_=mxGetField(f,0,"eval_pps"); \
    }

    // y=f(x)
    virtual void operator () (const X_Vector& x,Y_Vector& y) const { 
        X_Vector input[2]={eval_,x}; 
        mexCallMATLAB(1,&y,2,input,"feval");
    }

     // y=f'(x)dx 
     void p(
         const X_Vector& x,
         const X_Vector& dx,
         Y_Vector& y
     ) const { 
        X_Vector input[3]={eval_p_,x,dx}; 
        mexCallMATLAB(1,&y,3,input,"feval");
     }

     // z=f'(x)*dy
     void ps(
         const X_Vector& x,
         const Y_Vector& dy,
         X_Vector& z
     ) const {
        X_Vector input[3]={eval_ps_,x,dy}; 
        mexCallMATLAB(1,&z,3,input,"feval");
     }
     
     // z=(f''(x)dx)*dy
     void pps(
         const X_Vector& x,
         const X_Vector& dx,
         const Y_Vector& dy,
         X_Vector& z
     ) const { 
        X_Vector input[4]={eval_pps_,x,dx,dy}; 
        mexCallMATLAB(1,&z,4,input,"feval");
     }
};


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

    // Check that the structure contains all the functions
    for(int i=0;i<nops;i++)
        if(mxGetField(X,0,ops[i].c_str())==NULL)
            MatlabMessaging().error("Missing the " + ops[i] + 
                " function in the " + type + " " + name +".");
    
    // Check that all of the listed functions are really functions 
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
    peopt::ProblemClass::t problem_class = peopt::ProblemClass::Unconstrained;

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

// Check that we have a valid vector space.
void checkVS(const mxArray* VS) {
    // This implicitly checks that the vector space is valid
    peopt::ProblemClass::t problem_class_vs=get_problem_class_from_vs(VS);
}


// Check that we have a valid function bundle.
void checkFns(const mxArray* VS,const mxArray* fns) {
    // Infer problem classes for both the vector space and the functions
    peopt::ProblemClass::t problem_class_vs=get_problem_class_from_vs(VS);
    peopt::ProblemClass::t problem_class_fn=get_problem_class_from_fns(fns);

    // Make sure the problem class inferred by the vector space matches
    // that from the function class
    if(problem_class_vs != problem_class_fn) {
        std::stringstream ss;
        ss << "The problem class inferred from the vector space ("
            << peopt::ProblemClass::to_string(problem_class_vs) 
            << ") does not match\nthe problem class inferred from the "
            << "optimization functions ("
            << peopt::ProblemClass::to_string(problem_class_fn)
            << ").";

        MatlabMessaging().error(ss.str());
    }
}

// This sets the functions required for the vector spaces.
void setupVS(const mxArray* VS) { 
    // Get the problem class
    peopt::ProblemClass::t problem_class_vs=get_problem_class_from_vs(VS);

    // Setup the vector space
    switch(problem_class_vs) {
    case peopt::ProblemClass::Unconstrained:
        SET_VS(XX,mxGetField(VS,0,"X"));
        break;
    case peopt::ProblemClass::InequalityConstrained:
        SET_VS(XX,mxGetField(VS,0,"X"));
        SET_VS(ZZ,mxGetField(VS,0,"Z"));
        SET_EJA(ZZ,mxGetField(VS,0,"Z"));
        break;
    case peopt::ProblemClass::EqualityConstrained:
        SET_VS(XX,mxGetField(VS,0,"X"));
        SET_VS(YY,mxGetField(VS,0,"Y"));
        break;
    case peopt::ProblemClass::Constrained:
        SET_VS(XX,mxGetField(VS,0,"X"));
        SET_VS(YY,mxGetField(VS,0,"Y"));
        SET_VS(ZZ,mxGetField(VS,0,"Z"));
        SET_EJA(ZZ,mxGetField(VS,0,"Z"));
        break;
    }
}
// Check that the bundle of points aligns with the current vector space.
// We use this bundle of points for optimization 
void checkOptimizationPts(const mxArray* VS,const mxArray* pts) {
    // Get the problem class from the vector space
    peopt::ProblemClass::t problem_class=get_problem_class_from_vs(VS);
    
    // First, check that the bundle of points is a structure array        
    if(!mxIsStruct(pts))
        MatlabMessaging().error("The bundle of points must be a "
            "structure array.");

    // Next, insure that we have x, y, and z 
    if(mxGetField(pts,0,"x")==NULL)
        MatlabMessaging().error("The bundle of points must contain x.");
    if((problem_class == peopt::ProblemClass::EqualityConstrained 
        || problem_class == peopt::ProblemClass::Constrained) 
        && mxGetField(pts,0,"y")==NULL
    )
        MatlabMessaging().error("The bundle of points must contain y.");
    if((problem_class == peopt::ProblemClass::InequalityConstrained 
        || problem_class == peopt::ProblemClass::Constrained) 
        && mxGetField(pts,0,"z")==NULL
    )
        MatlabMessaging().error("The bundle of points must contain z.");
}

// Check that the bundle of points aligns with the current vector space.
// We use this bundle of points for the finite difference tests.
void checkDiagnosticPts(const mxArray* VS,const mxArray* pts) {
    // Get the problem class from the vector space
    peopt::ProblemClass::t problem_class=get_problem_class_from_vs(VS);

    // First, check that the bundle of points is a structure array        
    if(!mxIsStruct(pts))
        MatlabMessaging().error("The bundle of points must be a "
            "structure array.");

    // Next, insure that we have at least x, dx, and dxx 
    if(mxGetField(pts,0,"x")==NULL)
        MatlabMessaging().error("The bundle of points must contain x.");
    if(mxGetField(pts,0,"dx")==NULL)
        MatlabMessaging().error("The bundle of points must contain dx.");
    if(mxGetField(pts,0,"dxx")==NULL)
        MatlabMessaging().error("The bundle of points must contain dxx.");

    // Now, if we're equality constrained, check for dy 
    if(problem_class == peopt::ProblemClass::EqualityConstrained ||
       problem_class == peopt::ProblemClass::Constrained
    ) {
        if(mxGetField(pts,0,"dy")==NULL)
            MatlabMessaging().error("For equality constrained problems, "
                "the bundle of points must contain dy.");
    }

    // If we're inequality constrained, check for dz
    if(problem_class == peopt::ProblemClass::InequalityConstrained ||
       problem_class == peopt::ProblemClass::Constrained
    ) {
        if(mxGetField(pts,0,"dz")==NULL)
            MatlabMessaging().error("For inequality constrained problems, "
                "the bundle of points must contain dz.");
    }
}

// Check that our paramters are denoted by a string
void checkParams(const mxArray* params) {
    if(!mxIsChar(params))
        MatlabMessaging().error("The parameters must be a string that denotes "
            "a JSON parameter file.");
}

void mexFunction(
    int nOutput,mxArray* pOutput[],
    int nInput,const mxArray* pInput[]
) {
    // Create some type shortcuts
    typedef XX <double> X;
    typedef YY <double> Y;
    typedef ZZ <double> Z;
    
    // Check the number of arguments
    if(!(nInput==3 && nOutput==0) && !(nInput==4 && nOutput==1))
        MatlabMessaging().error("peopt usage:\n"
            "Optimization: x=peopt(VS,fns,pts,params)\n"
            "Diagnostics:  peopt(VS,fns,pts)");

    // Check the vector spaces
    checkVS(pInput[0]);

    // Check the bundle of functions
    checkFns(pInput[0],pInput[1]);

    // Determine if we're optimizing or doing diagnostics
    bool diagnostics = nInput==3;
     
    // Setup the vector space
    setupVS(pInput[0]);
   
    // Get the problem class
    peopt::ProblemClass::t problem_class=get_problem_class_from_vs(pInput[0]);

    // Create a bundle of functions
    peopt::Constrained <double,XX,YY,ZZ>::Functions::t fns;
    fns.f.reset(new MatlabScalarValuedFunction <double,XX>
        (mxGetField(pInput[1],0,"f")));
    switch(problem_class) {
    case peopt::ProblemClass::Unconstrained:
        break;
    case peopt::ProblemClass::EqualityConstrained:
        fns.g.reset(new MatlabVectorValuedFunction<double,XX,YY>
            (mxGetField(pInput[1],0,"g")));
        break;
    case peopt::ProblemClass::InequalityConstrained:
        fns.h.reset(new MatlabVectorValuedFunction<double,XX,ZZ>
            (mxGetField(pInput[1],0,"h")));
        break;
    case peopt::ProblemClass::Constrained:
        fns.g.reset(new MatlabVectorValuedFunction<double,XX,YY>
            (mxGetField(pInput[1],0,"g")));
        fns.h.reset(new MatlabVectorValuedFunction<double,XX,ZZ>
            (mxGetField(pInput[1],0,"h")));
        break;
    }

    // If we're doing diagnostics, diagnose
    if(diagnostics) {
        // Check that we have all the valid directions for the finite
        // difference tests.
        checkDiagnosticPts(pInput[0],pInput[2]);

        // Get the points for the objective finite difference tests
        X::Vector x=const_cast <X::Vector> (mxGetField(pInput[2],0,"x"));
        X::Vector dx=const_cast <X::Vector> (mxGetField(pInput[2],0,"dx"));
        X::Vector dxx=const_cast <X::Vector> (mxGetField(pInput[2],0,"dxx"));

        // Do a finite difference check and symmetric check on the objective
        MatlabMessaging().print("Diagnostics on the objective.");
        peopt::Diagnostics::gradientCheck<double,XX>
            (MatlabMessaging(),*(fns.f),x,dx);
        peopt::Diagnostics::hessianCheck <double,XX>
            (MatlabMessaging(),*(fns.f),x,dx);
        peopt::Diagnostics::hessianSymmetryCheck <double,XX>
            (MatlabMessaging(),*(fns.f),x,dx,dxx);

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == peopt::ProblemClass::EqualityConstrained ||
           problem_class == peopt::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the equality constraint.");
        
            // Get the points for the objective finite difference tests
            Y::Vector dy=const_cast <Y::Vector> (mxGetField(pInput[2],0,"dy"));

            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,XX,YY>
                (MatlabMessaging(),*(fns.g),x,dx,dy);
            peopt::Diagnostics::derivativeAdjointCheck <double,XX,YY>
                (MatlabMessaging(),*(fns.g),x,dx,dy);
            peopt::Diagnostics::secondDerivativeCheck <double,XX,YY>
                (MatlabMessaging(),*(fns.g),x,dx,dy);
        }

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == peopt::ProblemClass::InequalityConstrained ||
           problem_class == peopt::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the inequality constraint.");

            // Get the points for the objective finite difference tests
            Z::Vector dz=const_cast <Z::Vector> (mxGetField(pInput[2],0,"dz"));

            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,XX,ZZ>
                (MatlabMessaging(),*(fns.h),x,dx,dz);
            peopt::Diagnostics::derivativeAdjointCheck <double,XX,ZZ>
                (MatlabMessaging(),*(fns.h),x,dx,dz);
            peopt::Diagnostics::secondDerivativeCheck <double,XX,ZZ>
                (MatlabMessaging(),*(fns.h),x,dx,dz);
        }

    // Otherwise, let us solve an optimization problem
    } else {
        // Check that we have a valid parameter file
        checkParams(pInput[3]);
        mwSize buflen = mxGetNumberOfElements(pInput[3])+1;
        std::vector <char> params_(buflen);
        mxGetString(pInput[3],&(params_[0]),buflen);
        std::string params(params_.begin(),params_.end());
        
        // Check that we have all the valid directions for the finite
        // difference tests.
        checkOptimizationPts(pInput[0],pInput[2]);
        
        // Get the points for optimization 
        X::Vector x=const_cast <X::Vector> (mxGetField(pInput[2],0,"x"));

        // Do the optimization
        switch(problem_class) {
        case peopt::ProblemClass::Unconstrained: {
            peopt::Unconstrained <double,XX>::State::t state(x);
            peopt::json::Unconstrained <double,XX>
                ::read(MatlabMessaging(),params,state);
            peopt::Unconstrained<double,XX>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            X::init(state.x.back(),pOutput[0]);
            X::copy(state.x.back(),pOutput[0]);
            break;
        } case peopt::ProblemClass::InequalityConstrained: {
            Z::Vector z=const_cast <Z::Vector> (mxGetField(pInput[2],0,"z"));
            peopt::InequalityConstrained <double,XX,ZZ>::State::t state(x,z);
            peopt::json::InequalityConstrained <double,XX,ZZ>
                ::read(MatlabMessaging(),params,state);
            peopt::InequalityConstrained<double,XX,ZZ>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            X::init(state.x.back(),pOutput[0]);
            X::copy(state.x.back(),pOutput[0]);
            break;
        } case peopt::ProblemClass::EqualityConstrained: {
            MatlabMessaging().error(
                "Equality constrained optimization not currently implemented."); 
        } case peopt::ProblemClass::Constrained: {
            MatlabMessaging().error(
                "Fully constrained optimization not currently implemented."); 
        } }
    }
}
