#include <string>
#include <Python.h>
#include "peopt/peopt.h"
#include "peopt/json.h"

// A custom PyObject pointer that does proper clean-up on termination
struct PyObjectPtr { 
private:
    // Internal storage of the pointer
    PyObject* ptr;
    
public:
    // On default construction, we just create an empty pointer 
    PyObjectPtr() : ptr(NULL) {}

    // On construction, just initialize the auto_ptr
    PyObjectPtr(PyObject* ptr_) : ptr(ptr_) {}

    // Here, we copy both the pointer and increase the reference count
    // as long as the pointer is seemingly valid.
    PyObjectPtr& operator = (const PyObjectPtr& rhs) {
        ptr=rhs.get();
        if(ptr) Py_INCREF(ptr);
    }

    // For a reset, we decrement the pointer and then assign a new
    // value.
    void reset(PyObject* ptr_) {
        Py_XDECREF(ptr);
        ptr=ptr_;
    }

    // On a get, we simply return the pointer.
    PyObject* get() const {
        return ptr;
    }

    // On a release, we return the underlying pointer and then clear
    // the vector.  This will prevent a decrement later.
    PyObject* release() {
        PyObject* ptr_=ptr;
        ptr=NULL;
        return ptr_;
    }

    // On destruction, decrement the Python reference counter and do
    // not delete the pointer.
    ~PyObjectPtr() {
        Py_XDECREF(ptr);
        ptr=NULL;
    }

};

// A macro to alter the behavior of PyTuple_SetItem so that we don't
// have to hand increment the reference to the object since SetItem
// takes control of its arguments.
#define MyPyTuple_SetItem(p,pos,o) \
        Py_INCREF(o); \
        PyTuple_SetItem(p,pos,o);

struct PyVector : public PyObjectPtr {
private:
    // References to the algebra required in Python
    PyObjectPtr copy_; \
    PyObjectPtr scal_; \
    PyObjectPtr zero_; \
    PyObjectPtr axpy_; \
    PyObjectPtr innr_; \
    PyObjectPtr prod_; \
    PyObjectPtr id_; \
    PyObjectPtr linv_; \
    PyObjectPtr barr_; \
    PyObjectPtr srch_; \

public:
    // On basic initialization, just make sure that the internal storage is
    // NULL.
    PyVector() : PyObjectPtr(), copy_(NULL), scal_(NULL), zero_(NULL),
        axpy_(NULL), innr_(NULL), prod_(NULL), id_(NULL), linv_(NULL),
        barr_(NULL), srch_(NULL) {}

    // On a simple vector, initialize both the internal storage as well as
    // the basic linear algebra.
    PyVector(
        PyObject* vec,
        PyObject* copy__,
        PyObject* scal__,
        PyObject* zero__,
        PyObject* axpy__,
        PyObject* innr__
    ) : PyObjectPtr(vec), copy_(copy__), scal_(scal__), zero_(zero__),
        axpy_(axpy__), innr_(innr__), prod_(NULL), id_(NULL), linv_(NULL),
        barr_(NULL), srch_(NULL) {}

    // On a general vector, initialize both the internal storage as well as
    // all the linear algebra 
    PyVector(
        PyObject* vec,
        PyObject* copy__,
        PyObject* scal__,
        PyObject* zero__,
        PyObject* axpy__,
        PyObject* innr__,
        PyObject* prod__,
        PyObject* id__,
        PyObject* linv__,
        PyObject* barr__,
        PyObject* srch__
    ) : PyObjectPtr(vec), copy_(copy__), scal_(scal__), zero_(zero__),
        axpy_(axpy__), innr_(innr__), prod_(prod__), id_(id__), linv_(linv__),
        barr_(barr__), srch_(srch__) {}

    // Memory allocation and size setting 
    void init(const PyVector& vec) { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Grab the deepcopy function from the copy module 
        PyObjectPtr module(PyImport_ImportModule("copy")); 
        PyObjectPtr deepcopy(PyObject_GetAttrString(module.get(),"deepcopy")); 

        // Call deepcopy on vec and store the result internally
        PyObjectPtr args(PyTuple_New(1)); 
        MyPyTuple_SetItem(args.get(),0,vec.get()); 
        reset(PyObject_CallObject(deepcopy.get(),args.get())); 

        // Now, copy out all of the algebra functions
        copy_ = vec.copy_;
        scal_ = vec.scal_;
        zero_ = vec.zero_;
        axpy_ = vec.axpy_;
        innr_ = vec.innr_;
        prod_ = vec.prod_;
        id_ = vec.id_;
        linv_ = vec.linv_;
        barr_ = vec.barr_;
        srch_ = vec.srch_;

        // Release the global interpretter lock 
        PyGILState_Release(gstate);  
    } 
    
    // y <- x (Shallow.  No memory allocation.) 
    void copy(const PyVector& x) { 
        // Acquire the global interpreter lock
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the copy function on x and store internally 
        PyObjectPtr args(PyTuple_New(1)); 
        MyPyTuple_SetItem(args.get(),0,x.get()); 
        reset(PyObject_CallObject(copy_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    } 

    // x <- alpha * x
    void scal(const double& alpha) { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the scal function on alpha and the internal storage and store
        // in a temporary pointer.
        PyObjectPtr args(PyTuple_New(2)); 
        PyTuple_SetItem(args.get(),0,PyFloat_FromDouble(alpha)); 
        MyPyTuple_SetItem(args.get(),1,get()); 
        PyObject* alpha_x=PyObject_CallObject(scal_.get(),args.get()); 

        // Copy the result back into the internal storage 
        reset(alpha_x);

        /* Release the global interpretter lock */ 
        PyGILState_Release(gstate); 
    } 

    // x <- 0 
    void zero() { 
        // Acquire the global interpreter lock
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the zero function on this vector.  Store in z.
        PyObjectPtr args(PyTuple_New(1)); 
        MyPyTuple_SetItem(args.get(),0,get()); 
        PyObject* z=PyObject_CallObject(zero_.get(),args.get()); 

        // Copy the result back internally 
        reset(z);

        /* Release the global interpretter lock */ 
        PyGILState_Release(gstate); 
    } 

    // y <- alpha * x + y 
    void axpy(const double& alpha, const PyVector& x) { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the axpy function on alpha, x, and the internal storage.
        // Store in z.
        PyObjectPtr args(PyTuple_New(3)); 
        PyTuple_SetItem(args.get(),0,PyFloat_FromDouble(alpha)); 
        MyPyTuple_SetItem(args.get(),1,x.get()); 
        MyPyTuple_SetItem(args.get(),2,get()); 
        PyObject* z=PyObject_CallObject(axpy_.get(),args.get()); 

        // Copy the result back internally 
        reset(z); 

        /* Release the global interpretter lock */ 
        PyGILState_Release(gstate); 
    } 

    // innr <- <x,y> 
    double innr(const PyVector& x) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the innr function on x and the internal.  Store in zz. 
        PyObjectPtr args(PyTuple_New(2)); 
        MyPyTuple_SetItem(args.get(),0,x.get()); 
        MyPyTuple_SetItem(args.get(),1,get()); 
        PyObjectPtr zz(PyObject_CallObject(innr_.get(),args.get())); 

        // Grab a double of the result
        double z=PyFloat_AsDouble(zz.get()); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 

        // Return the result 
        return z; 
    } 

    // Jordan product, z <- x o y
    void prod(const PyVector& x, PyVector& z) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the prod function on x and the internal.  Store in z.
        PyObjectPtr args(PyTuple_New(2)); 
        MyPyTuple_SetItem(args.get(),0,x.get()); 
        MyPyTuple_SetItem(args.get(),1,get()); 
        z.reset(PyObject_CallObject(prod_.get(),args.get())); 

        /* Release the global interpretter lock */ 
        PyGILState_Release(gstate); 
    } 

    // Identity element, x <- e such that x o e = x 
    void id() { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the id function on the internal.  Store in z.
        PyObjectPtr args(PyTuple_New(1)); 
        MyPyTuple_SetItem(args.get(),0,get()); 
        PyObject* z=PyObject_CallObject(id_.get(),args.get()); 

        // Copy the result back into the internal. 
        reset(z);

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    } 

    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
    void linv(const PyVector& x, PyVector& z) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the linv function on x and the internal.  Store in z.
        PyObjectPtr args(PyTuple_New(2)); 
        MyPyTuple_SetItem(args.get(),0,x.get()); 
        MyPyTuple_SetItem(args.get(),1,get()); 
        z.reset(PyObject_CallObject(linv_.get(),args.get())); 

        // Release the global interpretter lock
        PyGILState_Release(gstate); 
    } 

    // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
    double barr() const { 
        // Acquire the global interpreter lock
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the barr function on the internal.  Store in zz.
        PyObjectPtr args(PyTuple_New(1)); 
        MyPyTuple_SetItem(args.get(),0,get()); 
        PyObjectPtr zz(PyObject_CallObject(barr_.get(),args.get())); 

        // Grab a double of the result 
        double z=PyFloat_AsDouble(zz.get()); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 

        // Return the result 
        return z; 
    } 

    // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} 
    // where y > 0.  If the argmax is infinity, then return Real(-1.). 
    double srch(const PyVector& x) const {  
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure(); 

        // Call the srch function on x and the internal.  Store in zz.
        PyObjectPtr args(PyTuple_New(2)); 
        MyPyTuple_SetItem(args.get(),0,x.get()); 
        MyPyTuple_SetItem(args.get(),1,get()); 
        PyObjectPtr zz(PyObject_CallObject(srch_.get(),args.get())); 

        // Grab a double of the result 
        double z=PyFloat_AsDouble(zz.get()); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 

        // Return the result 
        return z; 
    } 
};

template <typename Real=double> 
struct PythonVS { 
    // Setup the vector 
    typedef PyVector Vector; 

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
        y.prod(x,z);
    } 

    // Identity element, x <- e such that x o e = x 
    static void id(Vector& x) { 
        x.id();
    } 

    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
    static void linv(const Vector& x, const Vector& y, Vector& z) { 
        y.linv(x,z);
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
}; 

// A simple scalar valued function interface, f : X -> R
struct PythonScalarValuedFunction
    : peopt::ScalarValuedFunction <double,PythonVS>
{
private:
    // Python functions
    PyObjectPtr eval_;
    PyObjectPtr grad_;
    PyObjectPtr hessvec_;

    // Disallow the default constructor
    PythonScalarValuedFunction() : eval_(NULL), grad_(NULL), hessvec_(NULL) {}

public:
    // When we construct the function, we read in all the names from
    // the Python class.  Make sure that all of these functions are
    // valid prior to constructing this class.
    PythonScalarValuedFunction(PyObject* f)
        : eval_(NULL), grad_(NULL), hessvec_(NULL)
    {
        eval_.reset(PyObject_GetAttrString(f,"eval"));
        grad_.reset(PyObject_GetAttrString(f,"grad"));
        hessvec_.reset(PyObject_GetAttrString(f,"hessvec"));
    }

    // <- f(x) 
    double operator () (const PyVector& x) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the objective function on x 
        PyObjectPtr args(PyTuple_New(1));
        MyPyTuple_SetItem(args.get(),0,x.get());
        PyObjectPtr zz(PyObject_CallObject(eval_.get(),args.get())); 

        // Grab a double of the result and free the memory with the 
        // Python number.
        double z=PyFloat_AsDouble(zz.get());
       
        // Release the global interpretter lock 
        PyGILState_Release(gstate); 

        // Return the result
        return z;
    }

    // g = grad f(x) 
    void grad(const PyVector& x,PyVector& g) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the gradient function on x.  Return the answer in g
        PyObjectPtr args(PyTuple_New(1));
        MyPyTuple_SetItem(args.get(),0,x.get());
        g.reset(PyObject_CallObject(grad_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }

    // H_dx = hess f(x) dx 
    void hessvec(const PyVector& x,const PyVector& dx,PyVector& H_dx) const {
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the hessvec function on x and dx.  Return the answer in H_dx 
        PyObjectPtr args(PyTuple_New(2));
        MyPyTuple_SetItem(args.get(),0,x.get());
        MyPyTuple_SetItem(args.get(),1,dx.get());
        H_dx.reset(PyObject_CallObject(hessvec_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }
};

// A simple vector valued function interface, f : X -> Y
struct PythonVectorValuedFunction
    : public peopt::VectorValuedFunction<double,PythonVS,PythonVS> {
private:
    // Create some type shortcuts
    typedef PythonVS <>::Vector X_Vector; 
    typedef PythonVS <>::Vector Y_Vector; 

    // Python functions
    PyObjectPtr eval_;
    PyObjectPtr p_;
    PyObjectPtr ps_;
    PyObjectPtr pps_;

    // Prevent the default constructor from being called
    PythonVectorValuedFunction() : 
        eval_(NULL), p_(NULL), ps_(NULL), pps_(NULL) {}

public:
    // When we construct the function, we read in all the names from
    // the Python class.  Make sure that all of these functions are
    // valid prior to constructing this class.
    PythonVectorValuedFunction(PyObject* g) :
        eval_(NULL), p_(NULL), ps_(NULL), pps_(NULL)
    {
        eval_.reset(PyObject_GetAttrString(g,"eval"));
        p_.reset(PyObject_GetAttrString(g,"p"));
        ps_.reset(PyObject_GetAttrString(g,"ps"));
        pps_.reset(PyObject_GetAttrString(g,"pps"));
    }

    // y=f(x)
    void operator () (const X_Vector& x,Y_Vector& y) const {
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the evaluate function on x.  Return the answer in y.
        PyObjectPtr args(PyTuple_New(1));
        MyPyTuple_SetItem(args.get(),0,x.get());
        y.reset(PyObject_CallObject(eval_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }

    // y=f'(x)dx 
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Y_Vector& y
    ) const {
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the prime function on x and dx.  Return the answer in y.
        PyObjectPtr args(PyTuple_New(2));
        MyPyTuple_SetItem(args.get(),0,x.get());
        MyPyTuple_SetItem(args.get(),1,dx.get());
        y.reset(PyObject_CallObject(p_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }

    // z=f'(x)*dy
    virtual void ps(
        const X_Vector& x,
        const Y_Vector& dy,
        X_Vector& z
    ) const {
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the prime-adjoint function on x and dy.  Return the answer in z.
        PyObjectPtr args(PyTuple_New(2));
        MyPyTuple_SetItem(args.get(),0,x.get());
        MyPyTuple_SetItem(args.get(),1,dy.get());
        z.reset(PyObject_CallObject(ps_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }
     
    // z=(f''(x)dx)*dy
    virtual void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Y_Vector& dy,
        X_Vector& z
    ) const { 
        // Acquire the global interpreter lock 
        PyGILState_STATE gstate = PyGILState_Ensure();

        // Call the prime-adjoint function on x, dx, and dy.  Return the
        // answer in z.
        PyObjectPtr args(PyTuple_New(3));
        MyPyTuple_SetItem(args.get(),0,x.get());
        MyPyTuple_SetItem(args.get(),1,dx.get());
        MyPyTuple_SetItem(args.get(),2,dy.get());
        z.reset(PyObject_CallObject(pps_.get(),args.get())); 

        // Release the global interpretter lock 
        PyGILState_Release(gstate); 
    }
};

// A messaging utility that hooks directly into Matlab
struct PythonMessaging : public peopt::Messaging {
    // Prints a message
    void print(const std::string msg) const {
        PySys_WriteStdout((msg + '\n').c_str());
    }

    // Prints an error
    void error(const std::string msg) const {
        PySys_WriteStderr((msg + '\n').c_str());
        PyErr_SetString(PyExc_RuntimeError,msg.c_str());
        throw -1;
    }
};

// A finite difference test based on the vector space vs, functions fns, and
// points pts.
extern "C" PyObject* pypeopt(
    PyObject* opt_type_,
    PyObject* vs_,
    PyObject* fns_,
    PyObject* pts_,
    PyObject* fname_
) {
    // Acquire the global interpreter lock 
    PyGILState_STATE gstate = PyGILState_Ensure();
            
    try {
        // Try to read in a file name for the peopt parameters.  If it doesn't
        // exist, we assume that we're doing a finite difference test. 
        char* fname__ = PyString_AsString(fname_);
        std::string fname;
        if(fname__) fname.assign(fname__);
            
        // Determine the mode we're running in: 0=finite difference, 1=optimize 
        long mode = fname.size()==0 ? 0 : 1;

        // Determine the type of optimization
        // 0 = Unconstrained
        // 1 = Equality constrained
        // 2 = Inequality constrained
        // 3 = Constrained
        long opt_type = PyInt_AsLong(opt_type_);

        // Create the bundle of functions.
        peopt::Constrained <double,PythonVS,PythonVS,PythonVS>::Functions::t
            fns;

        // Add f, then see if we need to add g and h
        fns.f.reset(new PythonScalarValuedFunction (
            PyObject_GetAttrString(fns_,"f")));
        if(opt_type == 1 || opt_type == 3) 
            fns.g.reset(new PythonVectorValuedFunction (
                PyObject_GetAttrString(fns_,"g")));
        if(opt_type == 2 || opt_type == 3) 
            fns.h.reset(new PythonVectorValuedFunction(
                PyObject_GetAttrString(fns_,"h")));

        // Initialize points for x, y, and z when required.  In addition,
        // initialize points for dx, dxx, dy, and dz when doing the finite
        // difference tests
        PyObjectPtr X(PyObject_GetAttrString(vs_,"X"));
        std::auto_ptr <PyVector> x;
        std::auto_ptr <PyVector> dx;
        std::auto_ptr <PyVector> dxx;
        std::auto_ptr <PyVector> y;
        std::auto_ptr <PyVector> dy;
        std::auto_ptr <PyVector> z;
        std::auto_ptr <PyVector> dz;

        // Always initialize x
        x.reset(new PyVector(
            PyObject_GetAttrString(pts_,"x"),
            PyObject_GetAttrString(X.get(),"copy"),
            PyObject_GetAttrString(X.get(),"scal"),
            PyObject_GetAttrString(X.get(),"zero"),
            PyObject_GetAttrString(X.get(),"axpy"),
            PyObject_GetAttrString(X.get(),"innr")));

        // If we need to do a finite difference, initialize dx and dxx
        if(mode == 0) {
            dx.reset(new PyVector()); 
                dx->init(*x);
                dx->reset(PyObject_GetAttrString(pts_,"dx"));
            dxx.reset(new PyVector());
                dxx->init(*x);
                dxx->reset(PyObject_GetAttrString(pts_,"dxx"));
        }

        // Possibly initialize y and dy
        if(opt_type == 1 || opt_type == 3) {
            PyObjectPtr Y(PyObject_GetAttrString(vs_,"Y"));
            y.reset(new PyVector(
                PyObject_GetAttrString(pts_,"y"),
                PyObject_GetAttrString(Y.get(),"copy"),
                PyObject_GetAttrString(Y.get(),"scal"),
                PyObject_GetAttrString(Y.get(),"zero"),
                PyObject_GetAttrString(Y.get(),"axpy"),
                PyObject_GetAttrString(Y.get(),"innr")));
            if(mode == 0) 
                dy.reset(new PyVector()); 
                    dy->init(*y);
                    dy->reset(PyObject_GetAttrString(pts_,"dy"));
        }
        
        // Possibly initialize z and dz
        if(opt_type == 2 || opt_type == 3) {
            PyObjectPtr Z(PyObject_GetAttrString(vs_,"Z"));
            y.reset(new PyVector(
                PyObject_GetAttrString(pts_,"z"),
                PyObject_GetAttrString(Z.get(),"copy"),
                PyObject_GetAttrString(Z.get(),"scal"),
                PyObject_GetAttrString(Z.get(),"zero"),
                PyObject_GetAttrString(Z.get(),"axpy"),
                PyObject_GetAttrString(Z.get(),"innr"),
                PyObject_GetAttrString(Z.get(),"prod"),
                PyObject_GetAttrString(Z.get(),"id"),
                PyObject_GetAttrString(Z.get(),"linv"),
                PyObject_GetAttrString(Z.get(),"barr"),
                PyObject_GetAttrString(Z.get(),"srch")));
            if(mode == 0) 
                dz.reset(new PyVector()); 
                    dz->init(*z);
                    dz->reset(PyObject_GetAttrString(pts_,"dz"));
        }

        // Next, do the finite difference tests if required
        if(mode == 0) {
            // Do a finite difference check and symmetric check on the objective
            PythonMessaging().print("Diagnostics on the objective.");
            peopt::Diagnostics::gradientCheck<double,PythonVS>
                (PythonMessaging(),*(fns.f),*x,*dx);
            peopt::Diagnostics::hessianCheck <double,PythonVS>
                (PythonMessaging(),*(fns.f),*x,*dx);
            peopt::Diagnostics::hessianSymmetryCheck <double,PythonVS>
                (PythonMessaging(),*(fns.f),*x,*dx,*dxx);
        
            // Next, check if we have a equality constrained problem
            if(opt_type == 1 || opt_type == 3) {
                // Do some finite difference tests on the constraint
                peopt::Diagnostics::derivativeCheck <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.g),*x,*dx,*dy);
                peopt::Diagnostics::derivativeAdjointCheck
                    <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.g),*x,*dx,*dy);
                peopt::Diagnostics::secondDerivativeCheck
                    <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.g),*x,*dx,*dy);
            }

            // Finally, check if we have a inequality constrained problem
            if(opt_type == 2 || opt_type == 3) {
                // Do some finite difference tests on the constraint.
                peopt::Diagnostics::derivativeCheck
                    <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.h),*x,*dx,*dz);
                peopt::Diagnostics::derivativeAdjointCheck
                    <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.h),*x,*dx,*dz);
                peopt::Diagnostics::secondDerivativeCheck
                    <double,PythonVS,PythonVS>
                    (PythonMessaging(),*(fns.h),*x,*dx,*dz);
            }
        
            // Release the global interpretter lock 
            PyGILState_Release(gstate); 

            // Return none
            return Py_None;

        // Alternatively, optimize if required
        } else if(mode==1) {
            // Initialize memory for the solution
            PyObject* sol;

            // Optimize the unconstrained problem
            if(opt_type==0) {
                peopt::Unconstrained <double,PythonVS>::State::t state(*x);
                peopt::json::Unconstrained <double,PythonVS>::read(
                    PythonMessaging(),fname,state);
                peopt::Unconstrained <double,PythonVS>::Algorithms
                    ::getMin(PythonMessaging(),fns,state);
                sol = state.x.front().release();

            } else if(opt_type==1) {
                peopt::EqualityConstrained <double,PythonVS,PythonVS>::State::t
                    state(*x,*y);
                peopt::json::EqualityConstrained <double,PythonVS,PythonVS>
                    ::read(PythonMessaging(),fname,state);
                #if 0
                peopt::EqualityConstrained <double,PythonVS,PythonVS>
                    ::Algorithms::getMin(PythonMessaging(),fns,state);
                #endif
                sol = state.x.front().release();

            } else if(opt_type==2) {
                peopt::InequalityConstrained<double,PythonVS,PythonVS>::State::t
                    state(*x,*z);
                peopt::json::InequalityConstrained <double,PythonVS,PythonVS>
                    ::read(PythonMessaging(),fname,state);
                peopt::InequalityConstrained <double,PythonVS,PythonVS>
                    ::Algorithms::getMin(PythonMessaging(),fns,state);
                sol = state.x.front().release();

            } else if(opt_type==3) {
                peopt::Constrained<double,PythonVS,PythonVS,PythonVS>::State::t
                    state(*x,*y,*z);
                peopt::json::Constrained <double,PythonVS,PythonVS,PythonVS>
                    ::read(PythonMessaging(),fname,state);
                #if 0
                peopt::Constrained <double,PythonVS,PythonVS,PythonVS>
                    ::Algorithms::getMin(PythonMessaging(),fns,state);
                #endif
                sol = state.x.front().release();
            }

            // Release the global interpretter lock and return the solution
            PyGILState_Release(gstate); 
            return sol;
        }

    // In the case of an exception, the correct code should already be
    // set, so just exit.
    } catch (...) {
        // Release the global interpretter lock 
        PyGILState_Release(gstate); 

        // Return nothing
        return Py_None;
    }
}
