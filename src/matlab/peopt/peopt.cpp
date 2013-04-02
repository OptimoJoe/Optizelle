#include <string>
#include <sstream>
#include <memory>
#include "mex.h"
#include "peopt/peopt.h"
#include "peopt/json.h"

// Defines a simple Matlab autopointer
struct mxArrayPtr {
private:
    // Internal storage for the pointer
    mxArray* ptr;

public:
    // On the default construction, we create an empty pointer
    mxArrayPtr() : ptr(NULL) { }

    // On construction, just initialize the internal pointer
    mxArrayPtr(mxArray* ptr_) : ptr(ptr_) {}

    // For a reset, we free the memory and then assign the new vector
    void reset(mxArray* ptr_) {
        if(ptr) mxDestroyArray(ptr);
        ptr=ptr_;
    }

    // On a get, we return the pointer
    mxArray* get() const {
        return ptr;
    }

    // On a release, return the underlying pointer and then clear the
    // internal representation.
    mxArray* release() {
        mxArray* ptr_=ptr;
        ptr=NULL;
        return ptr_;
    }

    // On destruction, free the memory if it has been allocated
    ~mxArrayPtr() {
        if(ptr!=NULL) mxDestroyArray(ptr);
        ptr=NULL;
    }
};

// Handles a Matlab exception
void handleException(mxArray* err,std::string msg) {
    if(err) {
        // In the cass of an exception, grab the report
        mxArray* input[1]={err};
        mxArray* output[1];
        mexCallMATLAB(1,output,1,input,"getReport");

        // Turn the report into a string
        mwSize char_limit=256;
        char report_[char_limit];
        mxGetString(output[0],report_,char_limit);

        // The report has extra information that we don't want.  Hence,
        // we eliminate both the first line as well as the last two lines.
        // The first line is supposed to say what function this occured in,
        // but Matlab gets a little confused since we're doing mex trickery.
        // The last two lines will automatically be repeated by mexErrMsgTxt.
        std::string report=report_;
        size_t pos=report.find("\n");
        report = report.substr(pos+1);
        pos = report.rfind("\n");
        report = report.substr(0,pos);
        pos = report.rfind("\n");
        report = report.substr(0,pos);
        pos = report.rfind("\n");
        report = report.substr(0,pos);

        // Now, tack on our additional error message and then return control
        // to Matlab.
        std::stringstream ss;
        ss << msg << std::endl << std::endl << report;
        mexErrMsgTxt(ss.str().c_str());
    }
}

// Defines a vector for use in a vector space 
struct MatVector : public mxArrayPtr {
private:
    // References to the algebra required in Matlab
    mxArrayPtr copy_;
    mxArrayPtr scal_;
    mxArrayPtr zero_;
    mxArrayPtr axpy_;
    mxArrayPtr innr_;
    mxArrayPtr prod_;
    mxArrayPtr id_;
    mxArrayPtr linv_;
    mxArrayPtr barr_;
    mxArrayPtr srch_;
    mxArrayPtr symm_;

public:
    // On basic initialization, just make sure that the internal storage is
    // NULL.
    MatVector() : mxArrayPtr(), copy_(NULL), scal_(NULL), zero_(NULL),
        axpy_(NULL), innr_(NULL), prod_(NULL), id_(NULL), linv_(NULL),
        barr_(NULL), srch_(NULL), symm_(NULL) {}

    // On a simple vector, initialize both the internal storage as well as
    // the basic linear algebra.
    MatVector(
        mxArray* vec,
        mxArray* copy__,
        mxArray* scal__,
        mxArray* zero__,
        mxArray* axpy__,
        mxArray* innr__
    ) : prod_(NULL), id_(NULL), linv_(NULL), barr_(NULL), srch_(NULL),
        symm_(NULL)
    {
        // Compute a deep copy of the data in the supplied vector 
        reset(mxDuplicateArray(vec));

        // Now, copy out all of the algebra functions
        copy_.reset(mxDuplicateArray(copy__));
        scal_.reset(mxDuplicateArray(scal__));
        zero_.reset(mxDuplicateArray(zero__));
        axpy_.reset(mxDuplicateArray(axpy__));
        innr_.reset(mxDuplicateArray(innr__));
    }

    // On a general vector, initialize both the internal storage as well as
    // all the linear algebra 
    MatVector(
        mxArray* vec,
        mxArray* copy__,
        mxArray* scal__,
        mxArray* zero__,
        mxArray* axpy__,
        mxArray* innr__,
        mxArray* prod__,
        mxArray* id__,
        mxArray* linv__,
        mxArray* barr__,
        mxArray* srch__,
        mxArray* symm__
    ) { 
        // Compute a deep copy of the data in the supplied vector 
        reset(mxDuplicateArray(vec));

        // Now, copy out all of the algebra functions
        copy_.reset(mxDuplicateArray(copy__));
        scal_.reset(mxDuplicateArray(scal__));
        zero_.reset(mxDuplicateArray(zero__));
        axpy_.reset(mxDuplicateArray(axpy__));
        innr_.reset(mxDuplicateArray(innr__));
        prod_.reset(mxDuplicateArray(prod__));
        id_.reset(mxDuplicateArray(id__));
        linv_.reset(mxDuplicateArray(linv__));
        barr_.reset(mxDuplicateArray(barr__));
        srch_.reset(mxDuplicateArray(srch__));
        symm_.reset(mxDuplicateArray(symm__));
    }

    // Memory allocation and size setting 
    void init(const MatVector& x) {
        // Compute a deep copy of the data in x
        reset(mxDuplicateArray(x.get()));

        // Now, copy out all of the algebra functions
        copy_.reset(mxDuplicateArray(x.copy_.get()));
        scal_.reset(mxDuplicateArray(x.scal_.get()));
        zero_.reset(mxDuplicateArray(x.zero_.get()));
        axpy_.reset(mxDuplicateArray(x.axpy_.get()));
        innr_.reset(mxDuplicateArray(x.innr_.get()));
        if(x.prod_.get()) prod_.reset(mxDuplicateArray(x.prod_.get()));
        if(x.id_.get()) id_.reset(mxDuplicateArray(x.id_.get()));
        if(x.linv_.get()) linv_.reset(mxDuplicateArray(x.linv_.get()));
        if(x.barr_.get()) barr_.reset(mxDuplicateArray(x.barr_.get()));
        if(x.srch_.get()) srch_.reset(mxDuplicateArray(x.srch_.get()));
        if(x.symm_.get()) symm_.reset(mxDuplicateArray(x.symm_.get()));
    }

    // y <- x (Shallow.  No memory allocation.) 
    void copy(const mxArray* x){
        reset(mxDuplicateArray(x));
    }

    // y <- x (Shallow.  No memory allocation.) 
    void copy(const MatVector& x) {
        mxArray* input[2]={copy_.get(),x.get()};
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,"Error in the vector space copy function");
        reset(output[0]);
    }

    // x <- alpha * x 
    void scal(const double& alpha_) {
        // Create memory for the scalar
        mxArrayPtr alpha(mxCreateDoubleMatrix(1,1,mxREAL));
        mxGetPr(alpha.get())[0]=alpha_;

        // Create the inputs and outputs
        mxArray* input[3]={scal_.get(),alpha.get(),get()};
        mxArray* output[1];

        // Compute the scalar multiplication and store the result
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the vector space scal function");
        reset(output[0]);
    }

    // x <- 0 
    void zero() {
        // Create the inputs and outputs 
        mxArray* input[2]={zero_.get(),get()};
        mxArray* output[1];

        // Find the zero vector and store the result.
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,"Error in the vector space zero function");
        reset(output[0]);

    }

    // y <- alpha * x + y 
    void axpy(const double& alpha_, const MatVector& x) {
        // Create memory for the scalar
        mxArrayPtr alpha(mxCreateDoubleMatrix(1,1,mxREAL));
        mxGetPr(alpha.get())[0]=alpha_;

        // Create the inputs and outputs 
        mxArray* input[4]={axpy_.get(),alpha.get(),x.get(),get()};
        mxArray* output[1];

        // Compute the addition and store the result
        mxArray* err=mexCallMATLABWithTrap(1,output,4,input,"feval");
        handleException(err,"Error in the vector space axpy function");
        reset(output[0]);
    }

    // innr <- <x,y> 
    double innr(const MatVector& x) const {
        // Create the inputs and outputs 
        mxArray* input[3]={innr_.get(),x.get(),get()};
        mxArray* output[1];

        // Compute the inner product 
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the vector space innr function");
        mxArrayPtr output_(output[0]);

        // Get the result of the computation 
        double alpha=mxGetScalar(output_.get());

        // Return the result 
        return alpha;
    }

    // Jordan product, z <- x o y 
    void prod(const MatVector& x, const MatVector& y) {
        // Create the inputs and outputs 
        mxArray* input[3]={prod_.get(),x.get(),y.get()};
        mxArray* output[1];

        // Compute the product and store the result 
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the vector space prod function");
        reset(output[0]);
    }

    // Identity element, x <- e such that x o e = x 
    void id() {
        // Create the inputs and outputs 
        mxArray* input[2]={id_.get(),get()};
        mxArray* output[1];

        // Find the identity element and store the result
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,"Error in the vector space id function");
        reset(output[0]);
    }

    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
    void linv(const MatVector& x, const MatVector& y) {
        // Create the inputs and outputs 
        mxArray* input[3]={linv_.get(),x.get(),y.get()};
        mxArray* output[1];

        // Compute the product inverse and store the result 
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the vector space linv function");
        reset(output[0]);
    }

    // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
    double barr() const {
        // Create the inputs and outputs 
        mxArray* input[2]={barr_.get(),get()};
        mxArray* output[1];

        // Compute the barrier function 
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,"Error in the vector space barr function");
        mxArrayPtr output_(output[0]);

        // Get the result of the computation 
        double alpha=mxGetScalar(output_.get());

        // Return the result 
        return alpha;
    }

    // Line search, srch <- argmax {alphain Real >= 0 : alpha x + y >= 0} 
    // where y > 0.  If the argmax is infinity, then return Real(-1.). 
    double srch(const MatVector& x) const {
        // Create the inputs and outputs 
        mxArray* input[3]={srch_.get(),x.get(),get()};
        mxArray* output[1];

        // Compute the search function 
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the vector space srch function");
        mxArrayPtr output_(output[0]);
        
        // Get the result of the computation 
        double alpha=mxGetScalar(output_.get());

        // Return the result 
        return alpha;
    }

    // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
    // operator.
    void symm() {
        // Create the inputs and outputs 
        mxArray* input[2]={symm_.get(),get()};
        mxArray* output[1];

        // Find the symmetrization and store the result
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,"Error in the vector space symm function");
        reset(output[0]);
    }
};

        
template <typename Real=double> 
struct MatlabVS { 
    // Setup the vector 
    typedef MatVector Vector; 
        
    // Memory allocation and size setting 
    static void init(const Vector& x, Vector& y) {
        y.init(x);
    }       
                
    // y <- x (Shallow.  No memory allocation.) 
    static void copy(const Vector& x, Vector& y) {
        y.copy(x);
    } 
    
    // x <- alpha * x 
    static void scal(const Real& alpha, Vector& x) { 
        x.scal(alpha);
    }   
        
    // x <- 0 
    static void zero(Vector& x) { 
        x.zero();
    }
        
    // y <- alpha * x + y 
    static void axpy(const Real& alpha, const Vector& x, Vector& y) {
        y.axpy(alpha,x);
    }   
            
    // innr <- <x,y> 
    static Real innr(const Vector& x,const Vector& y) {
        return y.innr(x);
    } 
    
    // Jordan product, z <- x o y 
    static void prod(const Vector& x, const Vector& y, Vector& z) {
        z.prod(x,y);
    }   
        
    // Identity element, x <- e such that x o e = x 
    static void id(Vector& x) { 
        x.id();
    }   
        
    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
    static void linv(const Vector& x, const Vector& y, Vector& z) {
        z.linv(x,y);
    }
        
    // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
    static Real barr(const Vector& x) {
        return x.barr();
    }

    // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
    // where y > 0.  If the argmax is infinity, then return Real(-1.). 
    static Real srch(const Vector& x,const Vector& y) {
        return y.srch(x);
    }
        
    // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
    // operator.
    static void symm(Vector& x) {
        return x.symm();
    }
};

// A simple scalar valued function using Matlab provided functions 
struct MatlabScalarValuedFunction
    : public peopt::ScalarValuedFunction <double,MatlabVS>
{
private:
    // The Matlab functions
    mxArrayPtr eval_;
    mxArrayPtr grad_;
    mxArrayPtr hessvec_;

    // Disallow the default constructor
    MatlabScalarValuedFunction() : eval_(NULL), grad_(NULL), hessvec_(NULL) {}

public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabScalarValuedFunction(mxArray* f) {
        eval_.reset(mxDuplicateArray(mxGetField(f,0,"eval"))); 
        grad_.reset(mxDuplicateArray(mxGetField(f,0,"grad"))); 
        hessvec_.reset(mxDuplicateArray(mxGetField(f,0,"hessvec"))); 
    }

    // <- f(x) 
    double operator () (const MatVector& x) const {
        // Create the inputs and outputs 
        mxArray* input[2]={eval_.get(),x.get()}; 
        mxArray* output[1];

        // Evaluate the function 
        mexCallMATLAB(1,output,2,input,"feval"); 
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,
            "Error in the scalar valued function (objective) evaluation");
        mxArrayPtr output_(output[0]);

        // Get the result of the computation 
        double alpha=mxGetScalar(output_.get()); 

        // Return the result 
        return alpha; 
    }

    // g = grad f(x) 
    void grad(const MatVector& x,MatVector& g) const { 
        mxArray* input[2]={grad_.get(),x.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,
            "Error in the scalar valued function (objective) gradient");
        g.reset(output[0]);
    }

    // H_dx = hess f(x) dx 
    void hessvec(const MatVector& x,const MatVector& dx,MatVector& H_dx) const {
        mxArray* input[3]={hessvec_.get(),x.get(),dx.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,
            "Error in the scalar valued function (objective) Hessian-vector "
            "product");
        H_dx.reset(output[0]);
    }
};

// A vector valued function using Matlab provided functions 
struct MatlabVectorValuedFunction
    : public peopt::VectorValuedFunction <double,MatlabVS,MatlabVS> {
private:
    // Create some type shortcuts
    typedef MatlabVS <>::Vector X_Vector;
    typedef MatlabVS <>::Vector Y_Vector;

    // The Matlab functions
    mxArrayPtr eval_;
    mxArrayPtr p_;
    mxArrayPtr ps_;
    mxArrayPtr pps_;

    // Prevent the default constructor from being called
    MatlabVectorValuedFunction() :
        eval_(NULL), p_(NULL), ps_(NULL), pps_(NULL) {}

public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabVectorValuedFunction(mxArray* f) {
        eval_.reset(mxDuplicateArray(mxGetField(f,0,"eval")));
        p_.reset(mxDuplicateArray(mxGetField(f,0,"p"))); 
        ps_.reset(mxDuplicateArray(mxGetField(f,0,"ps")));
        pps_.reset(mxDuplicateArray(mxGetField(f,0,"pps"))); 
    }

    // y=f(x)
    void operator () (const X_Vector& x,Y_Vector& y) const { 
        mxArray* input[2]={eval_.get(),x.get()}; 
        mxArray* output[1];
        mexCallMATLAB(1,output,2,input,"feval");
        mxArray* err=mexCallMATLABWithTrap(1,output,2,input,"feval");
        handleException(err,
            "Error in the vector valued function (constraint) evaluation");
        y.reset(output[0]);
    }

     // y=f'(x)dx 
     void p(
         const X_Vector& x,
         const X_Vector& dx,
         Y_Vector& y
     ) const { 
        mxArray* input[3]={p_.get(),x.get(),dx.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,
            "Error in the vector valued function (constraint) derivative");
        y.reset(output[0]);
     }

     // z=f'(x)*dy
     void ps(
         const X_Vector& x,
         const Y_Vector& dy,
         X_Vector& z
     ) const {
        mxArray* input[3]={ps_.get(),x.get(),dy.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,
            "Error in the vector valued function (constraint) derivative "
            "adjoint");
        z.reset(output[0]);
     }
     
     // z=(f''(x)dx)*dy
     void pps(
         const X_Vector& x,
         const X_Vector& dx,
         const Y_Vector& dy,
         X_Vector& z
     ) const { 
        mxArray* input[4]={pps_.get(),x.get(),dx.get(),dy.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,4,input,"feval");
        handleException(err,
            "Error in the vector valued function (constraint) second "
            "derivative adjoint");
        z.reset(output[0]);
     }
};

// An operator using Matlab provided functions 
struct MatlabOperator : public peopt::Operator <double,MatlabVS,MatlabVS> {
private:
    // Create some type shortcuts
    typedef MatlabVS <>::Vector X_Vector;
    typedef MatlabVS <>::Vector Y_Vector;

    // The Matlab functions
    mxArrayPtr eval_;

    // The base point around where we evaluate the operator
    X_Vector& x; 
    
public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabOperator(mxArray* f,X_Vector& x_) : x(x_) {
        eval_.reset(mxDuplicateArray(mxGetField(f,0,"eval"))); 
    }

    // y=f(x)
    virtual void operator () (const X_Vector& dx,Y_Vector& y) const { 
        mxArray* input[3]={eval_.get(),x.get(),dx.get()}; 
        mxArray* output[1];
        mxArray* err=mexCallMATLABWithTrap(1,output,3,input,"feval");
        handleException(err,"Error in the operator");
        y.reset(output[0]);
    }
};

// A messaging utility that hooks directly into Matlab
struct MatlabMessaging : public peopt::Messaging {
    // Prints a message
    void print(const std::string msg) const {
        std::string msg_ = msg + "\n";
        mexPrintf(msg_.c_str());
        mexEvalString("drawnow");
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
    std::string ops[5] = { "copy", "scal", "zero", "axpy", "innr" };
    check_fns(X,ops,5,"vector space",name);
}

// Check that we have all the operations necessary for a Euclidean-Jordan
// algebra
void check_euclidean_jordan(const mxArray* X,std::string name) {
    // First, check that we have a valid vector space
    check_vector_space(X,name);

    // Next, we check we have a valid Euclidean-Jordan algebra
    std::string ops[6] = { "prod", "id", "linv", "barr", "srch", "symm" };
    check_fns(X,ops,6,"vector space",name);
}

// Check that we have all the operations necessary for a scalar valued function 
void check_scalar_valued_fn(const mxArray* X,std::string name) {
    std::string ops[3] = { "eval", "grad", "hessvec" };
    check_fns(X,ops,3,"scalar valued function",name);
}


// Check that we have all the operations necessary for a vector valued function 
void check_vector_valued_fn(const mxArray* X,std::string name) {
    std::string ops[4] = { "eval", "p", "ps", "pps" };
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
    get_problem_class_from_vs(VS);
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

// Check that the bundle of points aligns with the current vector space.
// We use this bundle of points for optimization 
void check_optimization_pts(const mxArray* VS,const mxArray* pts) {
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
void check_diagnostic_pts(const mxArray* VS,const mxArray* pts) {
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
     
    // Get the problem class
    peopt::ProblemClass::t problem_class=get_problem_class_from_vs(pInput[0]);

    // Create a bundle of functions
    peopt::Constrained <double,MatlabVS,MatlabVS,MatlabVS>::Functions::t fns;
    fns.f.reset(new MatlabScalarValuedFunction (mxGetField(pInput[1],0,"f")));
    switch(problem_class) {
    case peopt::ProblemClass::Unconstrained:
        break;
    case peopt::ProblemClass::EqualityConstrained:
        fns.g.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"g")));
        break;
    case peopt::ProblemClass::InequalityConstrained:
        fns.h.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"h")));
        break;
    case peopt::ProblemClass::Constrained:
        fns.g.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"g")));
        fns.h.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"h")));
        break;
    }

    // Initialize points for x, y, and z when required.  In addition,
    // initialize points for dx, dxx, dy, and dz when doing the finite
    // difference tests
    mxArray* X(mxGetField(pInput[0],0,"X"));
    std::auto_ptr <MatVector> x;
    std::auto_ptr <MatVector> dx;
    std::auto_ptr <MatVector> dxx;
    std::auto_ptr <MatVector> y;
    std::auto_ptr <MatVector> dy;
    std::auto_ptr <MatVector> z;
    std::auto_ptr <MatVector> dz;

    // Check that we have a valid set of optimization points (x,y,z)
    check_optimization_pts(pInput[0],pInput[2]);

    // Always initialize x
    x.reset(new MatVector(
        mxGetField(pInput[2],0,"x"),
        mxGetField(X,0,"copy"),
        mxGetField(X,0,"scal"),
        mxGetField(X,0,"zero"),
        mxGetField(X,0,"axpy"),
        mxGetField(X,0,"innr")));
    
    // If we need to do a finite difference, initialize dx and dxx
    if(diagnostics) {
        // Check that we have all the valid directions for the finite
        // difference tests.
        check_diagnostic_pts(pInput[0],pInput[2]);

        // Get the points for checking the objective
        dx.reset(new MatVector());
            dx->init(*x);
            dx->copy(mxGetField(pInput[2],0,"dx"));
        dxx.reset(new MatVector());
            dxx->init(*x);
            dxx->copy(mxGetField(pInput[2],0,"dx"));
    }
    
    // Possibly initialize y and dy
    if(problem_class == peopt::ProblemClass::EqualityConstrained ||
       problem_class == peopt::ProblemClass::Constrained
    ) {
        mxArray* Y(mxGetField(pInput[0],0,"Y"));
        y.reset(new MatVector(
            mxGetField(pInput[2],0,"y"),
            mxGetField(Y,0,"copy"),
            mxGetField(Y,0,"scal"),
            mxGetField(Y,0,"zero"),
            mxGetField(Y,0,"axpy"),
            mxGetField(Y,0,"innr")));
        if(diagnostics) {
            dy.reset(new MatVector());
                dy->init(*y);
                dy->copy(mxGetField(pInput[2],0,"dy"));
        }
    }

    // Possibly initialize z and dz
    if(problem_class == peopt::ProblemClass::InequalityConstrained ||
       problem_class == peopt::ProblemClass::Constrained
    ) {
        mxArray* Z(mxGetField(pInput[0],0,"Z"));
        z.reset(new MatVector(
            mxGetField(pInput[2],0,"z"),
            mxGetField(Z,0,"copy"),
            mxGetField(Z,0,"scal"),
            mxGetField(Z,0,"zero"),
            mxGetField(Z,0,"axpy"),
            mxGetField(Z,0,"innr"),
            mxGetField(Z,0,"prod"),
            mxGetField(Z,0,"id"),
            mxGetField(Z,0,"linv"),
            mxGetField(Z,0,"barr"),
            mxGetField(Z,0,"srch"),
            mxGetField(Z,0,"symm")));
        if(diagnostics) {
            dz.reset(new MatVector());
                dz->init(*z);
                dz->copy(mxGetField(pInput[2],0,"dz"));
        }
    }

    // If we're doing diagnostics, diagnose
    if(diagnostics) {
        // Do a finite difference check and symmetric check on the objective
        MatlabMessaging().print("Diagnostics on the objective.");
        peopt::Diagnostics::gradientCheck<double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx);
        peopt::Diagnostics::hessianCheck <double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx);
        peopt::Diagnostics::hessianSymmetryCheck <double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx,*dxx);

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == peopt::ProblemClass::EqualityConstrained ||
           problem_class == peopt::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the equality constraint.");
        
            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy);
            peopt::Diagnostics::derivativeAdjointCheck<double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy);
            peopt::Diagnostics::secondDerivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy);
        }

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == peopt::ProblemClass::InequalityConstrained ||
           problem_class == peopt::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the inequality constraint.");

            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz);
            peopt::Diagnostics::derivativeAdjointCheck<double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz);
            peopt::Diagnostics::secondDerivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz);
        }

    // Otherwise, let us solve an optimization problem
    } else {
        // Check that we have a valid parameter file
        checkParams(pInput[3]);
        mwSize buflen = mxGetNumberOfElements(pInput[3])+1;
        std::vector <char> params_(buflen);
        mxGetString(pInput[3],&(params_[0]),buflen);
        std::string params(params_.begin(),params_.end());

        // Create something to hold the solution
        const char* pnames[4]={"x","y","z","opt_stop"};
        pOutput[0]=mxCreateStructMatrix(1,1,4,pnames);
        
        // Do the optimization
        switch(problem_class) {
        case peopt::ProblemClass::Unconstrained: {
            // Allocate memory for the state
            peopt::Unconstrained <double,MatlabVS>::State::t state(*x);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x.back()));

            // Read the parameters and optimize
            peopt::json::Unconstrained <double,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            peopt::Unconstrained<double,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);

            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.back().release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                peopt::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case peopt::ProblemClass::InequalityConstrained: {
            // Allocate memory for the state
            peopt::InequalityConstrained <double,MatlabVS,MatlabVS>::State::t
                state(*x,*z);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x.back()));

            // Read the parameters and optimize
            peopt::json::InequalityConstrained <double,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            peopt::InequalityConstrained<double,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.back().release());
            mxSetField(pOutput[0],0,"z",state.z.back().release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                peopt::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case peopt::ProblemClass::EqualityConstrained: {
            // Allocate memory for the state
            peopt::EqualityConstrained <double,MatlabVS,MatlabVS>::State::t
                state(*x,*y);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x.back()));
            if(mxGetField(pInput[1],0,"PSchur_left")!=NULL)
                fns.PSchur_left.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_left"),state.x.back()));
            if(mxGetField(pInput[1],0,"PSchur_right")!=NULL)
                fns.PSchur_right.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_right"),state.x.back()));

            // Read the parameters and optimize
            peopt::json::EqualityConstrained <double,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            peopt::EqualityConstrained<double,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.back().release());
            mxSetField(pOutput[0],0,"y",state.y.back().release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                peopt::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case peopt::ProblemClass::Constrained: {
            // Allocate memory for the state
            peopt::Constrained <double,MatlabVS,MatlabVS,MatlabVS>::State::t
                state(*x,*y,*z);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x.back()));
            if(mxGetField(pInput[1],0,"PSchur_left")!=NULL)
                fns.PSchur_left.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_left"),state.x.back()));
            if(mxGetField(pInput[1],0,"PSchur_right")!=NULL)
                fns.PSchur_right.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_right"),state.x.back()));

            // Read the parameters and optimize
            peopt::json::Constrained <double,MatlabVS,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            peopt::Constrained<double,MatlabVS,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.back().release());
            mxSetField(pOutput[0],0,"y",state.y.back().release());
            mxSetField(pOutput[0],0,"y",state.z.back().release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                peopt::StoppingCondition::to_string(state.opt_stop).c_str()));
        } }
    }
}
