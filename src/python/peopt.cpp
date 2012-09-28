#include <Python.h>
#include "peopt.h"

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

    // On destruction, decrement the Python reference counter and do
    // not delete the pointer.
    ~PyObjectPtr() {
        Py_XDECREF(ptr);
    }
};

// A macro to alter the behavior of PyTuple_SetItem so that we don't
// have to hand increment the reference to the object since SetItem
// takes control of its arguments.
#define MyPyTuple_SetItem(p,pos,o) \
        Py_INCREF(o); \
        PyTuple_SetItem(p,pos,o);

// Define the Python vector space
#define CREATE_VS(name) \
template <typename Real> \
struct name { \
    /* Setup the vector */ \
    typedef PyObjectPtr Vector; \
\
    /* Store references to the Python algebra */ \
    static PyObjectPtr copy_; \
    static PyObjectPtr scal_; \
    static PyObjectPtr zero_; \
    static PyObjectPtr axpy_; \
    static PyObjectPtr innr_; \
    static PyObjectPtr prod_; \
    static PyObjectPtr id_; \
    static PyObjectPtr linv_; \
    static PyObjectPtr barr_; \
    static PyObjectPtr srch_; \
\
    /* Memory allocation and size setting */ \
    static void init(const Vector& x, Vector& y) { \
        /* Acquire the global interpreter lock */  \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Grab the deepcopy function from the copy module */ \
        PyObjectPtr module(PyImport_ImportModule("copy")); \
        PyObjectPtr deepcopy(PyObject_GetAttrString(module.get(),"deepcopy")); \
\
        /* Call deepcopy on x and store the result in y */ \
        PyObjectPtr args(PyTuple_New(1)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        y.reset(PyObject_CallObject(deepcopy.get(),args.get())); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate);  \
    } \
\
    /* y <- x (Shallow.  No memory allocation.) */ \
    static void copy(const Vector& x, Vector& y) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the copy function on x */ \
        PyObjectPtr args(PyTuple_New(1)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        y.reset(PyObject_CallObject(copy_.get(),args.get())); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* x <- alpha * x */ \
    static void scal(const Real& alpha, Vector& x) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the scal function on alpha and x */ \
        PyObjectPtr args(PyTuple_New(2)); \
        PyTuple_SetItem(args.get(),0,PyFloat_FromDouble(alpha)); \
        MyPyTuple_SetItem(args.get(),1,x.get()); \
        PyObjectPtr alpha_x(PyObject_CallObject(scal_.get(),args.get())); \
\
        /* Copy the result back into x and free the memory associated */ \
        /* with the temporary. */ \
        copy(alpha_x,x); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* x <- 0 */ \
    static void zero(Vector& x) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the zero function on x */ \
        PyObjectPtr args(PyTuple_New(1)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        PyObjectPtr z(PyObject_CallObject(zero_.get(),args.get())); \
\
        /* Copy the result back into x and free the memory associated */ \
        /* with the temporary. */ \
        copy(z,x); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* y <- alpha * x + y */ \
    static void axpy(const Real& alpha, const Vector& x, Vector& y) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the axpy function on alpha, x, and y */ \
        PyObjectPtr args(PyTuple_New(3)); \
        PyTuple_SetItem(args.get(),0,PyFloat_FromDouble(alpha)); \
        MyPyTuple_SetItem(args.get(),1,x.get()); \
        MyPyTuple_SetItem(args.get(),2,y.get()); \
        PyObjectPtr z(PyObject_CallObject(axpy_.get(),args.get())); \
\
        /* Copy the result back into y and free the memory associated */ \
        /* with the temporary */ \
        copy(z,y); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* innr <- <x,y> */ \
    static Real innr(const Vector& x,const Vector& y) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the innr function on x and y */ \
        PyObjectPtr args(PyTuple_New(2)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        MyPyTuple_SetItem(args.get(),1,y.get()); \
        PyObjectPtr zz(PyObject_CallObject(innr_.get(),args.get())); \
\
        /* Grab a double of the result and free the memory with the */ \
        /* Python number. */ \
        double z=PyFloat_AsDouble(zz.get()); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
\
        /* Return the result */ \
        return z; \
    } \
\
    /* Jordan product, z <- x o y */ \
    static void prod(const Vector& x, const Vector& y, Vector& z) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the prod function on x and y */ \
        PyObjectPtr args(PyTuple_New(2)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        MyPyTuple_SetItem(args.get(),1,y.get()); \
        z.reset(PyObject_CallObject(prod_.get(),args.get())); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* Identity element, x <- e such that x o e = x */  \
    static void id(Vector& x) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the id function on x */ \
        PyObjectPtr args(PyTuple_New(1)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        PyObjectPtr z(PyObject_CallObject(id_.get(),args.get())); \
\
        /* Copy the result back into x and free the memory associated */ \
        /* with the temporary. */ \
        copy(z,x); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y */ \
    static void linv(const Vector& x, const Vector& y, Vector& z) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the linv function on x and y */ \
        PyObjectPtr args(PyTuple_New(2)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        MyPyTuple_SetItem(args.get(),1,y.get()); \
        z.reset(PyObject_CallObject(linv_.get(),args.get())); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
    } \
\
    /* Barrier function, barr <- barr(x) where x o grad barr(x) = e */  \
    static Real barr(const Vector& x) { \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the barr function on x */ \
        PyObjectPtr args(PyTuple_New(1)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        PyObjectPtr zz(PyObject_CallObject(barr_.get(),args.get())); \
\
        /* Grab a double of the result and free the memory with the */ \
        /* Python number. */ \
        double z=PyFloat_AsDouble(zz.get()); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
\
        /* Return the result */ \
        return z; \
    } \
\
    /* Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >= 0} */ \
    /* where y > 0.  If the argmax is infinity, then return Real(-1.). */ \
    static Real srch(const Vector& x,const Vector& y) {  \
        /* Acquire the global interpreter lock */ \
        PyGILState_STATE gstate = PyGILState_Ensure(); \
\
        /* Call the srch function on x and y*/ \
        PyObjectPtr args(PyTuple_New(2)); \
        MyPyTuple_SetItem(args.get(),0,x.get()); \
        MyPyTuple_SetItem(args.get(),1,y.get()); \
        PyObjectPtr zz(PyObject_CallObject(srch_.get(),args.get())); \
\
        /* Grab a double of the result and free the memory with the */ \
        /* Python number. */ \
        double z=PyFloat_AsDouble(zz.get()); \
\
        /* Release the global interpretter lock */ \
        PyGILState_Release(gstate); \
\
        /* Return the result */ \
        return z; \
    } \
}; \
\
/* Initialize the starting pointer values for each of these classes */ \
template <typename Real> PyObjectPtr name <Real>::copy_(NULL); \
template <typename Real> PyObjectPtr name <Real>::scal_(NULL); \
template <typename Real> PyObjectPtr name <Real>::zero_(NULL); \
template <typename Real> PyObjectPtr name <Real>::axpy_(NULL); \
template <typename Real> PyObjectPtr name <Real>::innr_(NULL); \
template <typename Real> PyObjectPtr name <Real>::prod_(NULL); \
template <typename Real> PyObjectPtr name <Real>::id_(NULL); \
template <typename Real> PyObjectPtr name <Real>::linv_(NULL); \
template <typename Real> PyObjectPtr name <Real>::barr_(NULL); \
template <typename Real> PyObjectPtr name <Real>::srch_(NULL); 

// This is a shortcut to setup a vector space
#define SET_VS(name,obj,str,from) \
    PyObjectPtr obj(PyObject_GetAttrString(from,str)); \
    name <double>::copy_.reset(PyObject_GetAttrString(obj.get(),"copy")); \
    name <double>::scal_.reset(PyObject_GetAttrString(obj.get(),"scal")); \
    name <double>::zero_.reset(PyObject_GetAttrString(obj.get(),"zero")); \
    name <double>::axpy_.reset(PyObject_GetAttrString(obj.get(),"axpy")); \
    name <double>::innr_.reset(PyObject_GetAttrString(obj.get(),"innr"));

// This is a shortcut to setup a Euclidean-Jordan algebra
#define SET_EJA(name,obj,str,from) \
    PyObjectPtr obj(PyObject_GetAttrString(from,str)); \
    name <double>::prod_.reset(PyObject_GetAttrString(obj.get(),"prod")); \
    name <double>::id_.reset(PyObject_GetAttrString(obj.get(),"id")); \
    name <double>::linv_.reset(PyObject_GetAttrString(obj.get(),"linv")); \
    name <double>::barr_.reset(PyObject_GetAttrString(obj.get(),"barr")); \
    name <double>::srch_.reset(PyObject_GetAttrString(obj.get(),"srch"));

// Creat the vector spaces
CREATE_VS(XX);
CREATE_VS(YY);
CREATE_VS(ZZ);

// A simple scalar valued function interface, f : X -> R
template <typename Real,template <typename> class XX>
struct PythonScalarValuedFunction : peopt::ScalarValuedFunction <double,XX> {
private:
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector Vector;

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
    Real operator () (const Vector& x) const { 
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
    void grad(const Vector& x,Vector& g) const { 
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
    void hessvec(const Vector& x,const Vector& dx,Vector& H_dx) const {
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
template <
    typename Real,
    template <typename> class XX,
    template <typename> class YY 
>  
struct PythonVectorValuedFunction
    : public peopt::VectorValuedFunction<Real,XX,YY> {
private:
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector; 
    typedef YY <Real> Y;
    typedef typename Y::Vector Y_Vector; 

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
extern "C" void finite_difference_test(
    PyObject* opt_type_,
    PyObject* vs_,
    PyObject* fns_,
    PyObject* pts_
) {
    // Acquire the global interpreter lock 
    PyGILState_STATE gstate = PyGILState_Ensure();

    try {
        // Determine the type of optimization
        long opt_type = PyInt_AsLong(opt_type_);

        // Create the vector space X
        SET_VS(XX,X,"X",vs_);
        
        // Create the bundle of functions and add f
        peopt::Constrained <double,XX,YY,ZZ>::Functions::t fns;
        fns.f.reset(new PythonScalarValuedFunction <double,XX> (
            PyObject_GetAttrString(fns_,"f")));

        // Create the points for the finite difference test for the objective
        typedef XX <double>::Vector X_Vector;
        X_Vector x(PyObject_GetAttrString(pts_,"x"));
        X_Vector dx(PyObject_GetAttrString(pts_,"dx"));
        X_Vector dxx(PyObject_GetAttrString(pts_,"dxx"));

        // Do a finite difference check and symmetric check on the objective
        PythonMessaging().print("Diagnostics on the objective.");
        peopt::Diagnostics::gradientCheck<double,XX>
            (PythonMessaging(),*(fns.f),x,dx);
        peopt::Diagnostics::hessianCheck <double,XX>
            (PythonMessaging(),*(fns.f),x,dx);
        peopt::Diagnostics::hessianSymmetryCheck <double,XX>
            (PythonMessaging(),*(fns.f),x,dx,dxx);

        // Next, check if we have a equality constrained problem
        if(opt_type == 1 || opt_type == 3) {
            // Create the vector space Y
            SET_VS(YY,Y,"Y",vs_);
        
            // Add g to the bundle of functions
            fns.g.reset(new PythonVectorValuedFunction <double,XX,YY> (
                PyObject_GetAttrString(fns_,"g")));
        
            // Create the points for the finite difference test for the
            // constraint.
            typedef YY <double>::Vector Y_Vector;
            Y_Vector y(PyObject_GetAttrString(pts_,"y"));
            Y_Vector dy(PyObject_GetAttrString(pts_,"dy"));

            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,XX,YY>
                (PythonMessaging(),*(fns.g),x,dx,dy);
            peopt::Diagnostics::derivativeAdjointCheck <double,XX,YY>
                (PythonMessaging(),*(fns.g),x,dx,dy);
            peopt::Diagnostics::secondDerivativeCheck <double,XX,YY>
                (PythonMessaging(),*(fns.g),x,dx,dy);
        }

        // Finally, check if we have an inequality constrained problem
        if(opt_type == 2 || opt_type == 3) {
            // Create the vector space Z
            SET_VS(ZZ,Z,"Z",vs_);
            SET_EJA(ZZ,Z0,"Z",vs_);
        
            // Add g to the bundle of functions
            fns.h.reset(new PythonVectorValuedFunction <double,XX,ZZ> (
                PyObject_GetAttrString(fns_,"h")));
        
            // Create the points for the finite difference test for the
            // constraint.
            typedef ZZ <double>::Vector Z_Vector;
            Z_Vector z(PyObject_GetAttrString(pts_,"z"));
            Z_Vector dz(PyObject_GetAttrString(pts_,"dz"));

            // Do some finite difference tests on the constraint
            peopt::Diagnostics::derivativeCheck <double,XX,ZZ>
                (PythonMessaging(),*(fns.h),x,dx,dz);
            peopt::Diagnostics::derivativeAdjointCheck <double,XX,ZZ>
                (PythonMessaging(),*(fns.h),x,dx,dz);
            peopt::Diagnostics::secondDerivativeCheck <double,XX,ZZ>
                (PythonMessaging(),*(fns.h),x,dx,dz);
        }
    } catch (...) {}
        
    // Release the global interpretter lock 
    PyGILState_Release(gstate); 
}
