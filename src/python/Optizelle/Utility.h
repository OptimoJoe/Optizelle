#pragma once

#include <Python.h>
#include <optizelle/exception.h>
#include <optizelle/optizelle.h>
#include <optizelle/json.h>

// Alright, integrating C++ with Python is fraught with issues, but one that
// affects us specifically is const correctness.  There's not really a good
// way to handle this since Python doesn't have a concept of constant elements.
// In our code, we will consider something constant if we do not change its
// underlying value.  Now, we may attach to it, which increments its reference
// count, but we will consider this to be a constant operation.

namespace Optizelle {
    // Predeclare out pointer types
    namespace Python {
        struct PyObjectPtr;
    }

    namespace OptimizationStop {
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & opt_stop);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace TruncatedStop {
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & trunc_stop);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace AlgorithmClass {
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & algorithm_class);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace Operators{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & op);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace LineSearchDirection{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & dir);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace LineSearchKind{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & kind);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace OptimizationLocation{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & loc);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace FunctionDiagnostics{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & diag);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace VectorSpaceDiagnostics{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & diag);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace DiagnosticScheme{
        // Converts t to a Python enumerated type
        Python::PyObjectPtr toPython(t const & dscheme);

        // Converts a Python enumerated type to t
        t fromPython(Python::PyObjectPtr const & member);
    }

    namespace Python {
        // Exception for when a function in Python throws an error
        namespace Exception {
            struct t : public std::runtime_error {
                using std::runtime_error::runtime_error;
            };
        }

        // Manage the memory and reference counts for raw PyObject pointers
        struct PyObjectPtr {
        protected:
            // Internal storage of the pointer
            PyObject * ptr;

        public:
            // Determine whether we're accepting a borrowed reference or a new
            // reference
            enum Mode : std::size_t {
                Borrowed,       // Borrowed reference
                New             // New reference
            };

            // Grab the pointer
            PyObjectPtr(
                PyObject * const ptr_,
                PyObjectPtr::Mode const & mode = PyObjectPtr::New
            );

            // Copy semantics
            PyObjectPtr(PyObjectPtr const & p);
            PyObjectPtr & operator = (PyObjectPtr& p);

            // Move semantics
            PyObjectPtr(PyObjectPtr && p);
            PyObjectPtr & operator = (PyObjectPtr && p);

            // On a get, we simply return the pointer.
            PyObject * get() const;

            // On destruction, decrement the reference count
            ~PyObjectPtr();
        };

        // Wrappers to the Python C-API.  The goal here is to wrap these
        // with our memory management scheme.
        namespace capi {
            PyObjectPtr PyImport_ImportModule(const char *name);

            PyObjectPtr PyUnicode_FromString(const char *v);
            std::string PyUnicode_AsUTF8(PyObjectPtr const & string);

            PyObjectPtr PyInt_FromNatural(Natural const & ival);
            Natural PyInt_AsNatural(PyObjectPtr const & io);

            PyObjectPtr PyFloat_FromDouble(double v);
            double PyFloat_AsDouble(PyObjectPtr const & pyfloat);

            PyObjectPtr PyObject_GetAttrString(
                PyObjectPtr const & o,
                const char *attr_name);
            void PyObject_SetAttrString(
                PyObjectPtr & o,
                const char * attr_name,
                PyObjectPtr const & v);
            PyObjectPtr PyObject_CallObject1(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                std::string const & errmsg);
            PyObjectPtr PyObject_CallObject2(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                std::string const & errmsg);
            PyObjectPtr PyObject_CallObject3(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                PyObjectPtr const & arg3,
                std::string const & errmsg);
            PyObjectPtr PyObject_CallObject4(
                PyObjectPtr const & fn,
                PyObjectPtr const & arg1,
                PyObjectPtr const & arg2,
                PyObjectPtr const & arg3,
                PyObjectPtr const & arg4,
                std::string const & errmsg);

            PyObjectPtr PyTuple_New(Py_ssize_t const & len);
            void PyTuple_SetItem(
                PyObjectPtr const & p,
                Py_ssize_t const & pos,
                PyObjectPtr const & o);
            PyObjectPtr PyTuple_GetItem(
                PyObjectPtr const & p,
                Py_ssize_t const & pos);
            PyObjectPtr PyTuple_Pack_2(
                PyObjectPtr const & item1,
                PyObjectPtr const & item2);

            PyObjectPtr PyList_New(Py_ssize_t const & len);
            void PyList_Append(PyObjectPtr & list, PyObjectPtr const & item);
            Natural PyList_Size(PyObjectPtr const & list);
            PyObjectPtr PyList_GetItem(
                PyObjectPtr const & list,
                Py_ssize_t const & index);

            // Calls the Optizelle exception with a string
            void PyErr_SetString_Optizelle(std::string const & msg);

            // Deep copy of a Python object and return the result
            PyObjectPtr deepcopy(PyObjectPtr const & in);

            // Converts an Optizelle enumerated type to a PyObject
            PyObjectPtr enumToPyObject(
                std::string const & type,
                std::string const & member
            );

            // Converts an Optizelle enumerated type to a Natural based on
            // the scheme in the Python enumerated type
            Natural enumToNatural(
                std::string const & type,
                std::string const & member
            );
        }

        // A messaging utility that hooks directly into Python
        namespace Messaging {
            Optizelle::Messaging::t python(PyObjectPtr const & print);
        }

        // This class merges the vector space with a vector into a singular
        // object.  We require this structure since Optizelle requires the
        // vector space to be static.  Since the user is passing us a vector
        // space dynamically, we merge the vector space functions with the
        // vectors and then statically define the vector space to call these
        // functions.
        struct Vector {
            // Vector space
            PyObjectPtr vs;

            // Vector data
            PyObjectPtr data;

            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Vector)

            // Create a vector with the appropriate vector space
            Vector(PyObjectPtr const & vs_,PyObjectPtr const & vec_);

            // Move semantics
            Vector(Vector && vec) = default;
            Vector & operator = (Vector && vec) = default;

            // Memory allocation and size setting
            Vector init() const;

            // y <- x (Shallow.  No memory allocation.)  Internal is y.
            void copy(Vector const & x);

            // x <- alpha * x.  Internal is x.
            void scal(double const & alpha);

            // x <- 0.  Internal is x.
            void zero();

            // y <- alpha * x + y.  Internal is y.
            void axpy(double const & alpha, Vector const & x);

            // innr <- <x,y>.  Internal is y.
            double innr(Vector const & x) const;

            // x <- random.  Internal is x.
            void rand();

            // Jordan product, z <- x o y.  Internal is z.
            void prod(Vector const & x, Vector const & y);

            // Identity element, x <- e such that x o e = x.  Internal is x.
            void id();

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y .
            // Internal is z.
            void linv(Vector const & x, Vector const & y);

            // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
            // Internal is x.
            double barr() const;

            // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >=
            // 0} where y > 0.  Internal is y.
            double srch(Vector const & x) const;

            // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
            // operator.  Internal is x.
            void symm();

            // Converts (copies) a value into Python.  This assumes memory
            // has been allocated both in the vector as well as Python.
            void toPython(PyObjectPtr const & data) const;

            // Converts (copies) a value from Python.  This assumes memory
            // has been allocated both in the vector as well as Python.
            void fromPython(PyObjectPtr const & data);
        };

        // Python state
        template <typename ProblemClass>
        struct State {
            // Data
            PyObjectPtr data;

            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(State)

            // Add move semantics
            State(State &&) = default;
            State & operator = (State &&) = default;

            // On construction, we just grab the pointer to the state object
            State(
                PyObjectPtr const & data_
            ) : data(data_) {}

            // Convert a C++ state to a Python state
            void toPython(typename ProblemClass::State::t const & state);

            // Convert a Python state to C++
            void fromPython(typename ProblemClass::State::t & state);
        };

        // Python bundle of functions
        template <typename ProblemClass>
        struct Functions {
        private:
            // Keep some states lying around so that we can communicate this
            // to our operator.
            State <ProblemClass> & pystate;
            typename ProblemClass::State::t const & state;

        public:
            // Data
            PyObjectPtr data;

            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Functions)

            // Add move semantics
            Functions(Functions &&) = default;
            Functions & operator = (Functions &&) = default;

            // On construction, we just grab the pointer to the bundle object
            Functions(
                State <ProblemClass> & pystate_,
                typename ProblemClass::State::t const & state_,
                PyObjectPtr const & data_
            ) :
                pystate(pystate_),
                state(state_),
                data(data_)
            {}

            // Convert a Python bundle to C++
            void fromPython(typename ProblemClass::Functions::t & fns);
        };

        // The state manipulator for Python
        template <typename ProblemClass>
        struct StateManipulator :
            public Optizelle::StateManipulator <ProblemClass>
        {
        private:
            // Keep a copy of a Python state lying around so that we can
            // use it to pass information back and forth to the Python
            // StateManipulator
            State <ProblemClass> & pystate;

            // Similarly, we the bundle of functions lying around
            Functions <ProblemClass> & pyfns;

            // Data
            PyObjectPtr data;

        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(StateManipulator)

            // Add move semantics
            StateManipulator(StateManipulator && smanip) :
                pystate(smanip.pystate),
                pyfns(smanip.pyfns),
                data(std::move(smanip.data))
            {}

            // We need the Python state manipulator, a copy of a Python state
            // to pass information, and a copy of the Python functions.
            StateManipulator(
                State <ProblemClass> & pystate_,
                Functions <ProblemClass> & pyfns_,
                PyObjectPtr const & data_
            ) :
                pystate(pystate_),
                pyfns(pyfns_),
                data(data_)
            {}

            // Application
            void eval(
                typename ProblemClass::Functions::t const & fns,
                typename ProblemClass::State::t & state,
                OptimizationLocation::t const & loc_
            ) const {
                // Convert the C++ state to a Python state
                pystate.toPython(state);

                // Convert the lcoation to Python
                auto loc = OptimizationLocation::toPython(loc_);

                // Call the Python state manipulator
                auto eval = capi::PyObject_GetAttrString(data,"eval");
                capi::PyObject_CallObject3(
                    eval,
                    pyfns.data,
                    pystate.data,
                    loc,
                    __LOC__
                        + ", evaluation of the StateManipulator object failed");

                // Convert the Python state to the C++ state
                pystate.fromPython(state);
            }
        };

        // Vector space that works with Python objects
        template <typename Real=double>
        struct PythonVS {
            // Prevent constructors
            NO_CONSTRUCTORS(PythonVS)

            // Setup the vector
            typedef Optizelle::Python::Vector Vector;

            // Memory allocation and size setting
            static Vector init(Vector const & x) {
                return x.init();
            }

            // y <- x (Shallow.  No memory allocation.)
            static void copy(Vector const & x, Vector & y) {
                y.copy(x);
            }

            // x <- alpha * x
            static void scal(Real const & alpha, Vector & x) {
                x.scal(alpha);
            }

            // x <- 0
            static void zero(Vector & x) {
                x.zero();
            }

            // y <- alpha * x + y
            static void axpy(Real const & alpha, Vector const & x, Vector & y) {
                y.axpy(alpha,x);
            }

            // innr <- <x,y>
            static Real innr(Vector const & x,Vector const & y) {
                return y.innr(x);
            }

            // x <- random
            static void rand(Vector & x){
                x.rand();
            }

            // Jordan product, z <- x o y
            static void prod(Vector const & x, Vector const & y, Vector & z) {
                z.prod(x,y);
            }

            // Identity element, x <- e such that x o e = x
            static void id(Vector & x) {
                x.id();
            }

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
            static void linv(Vector const & x, Vector const & y, Vector & z) {
                z.linv(x,y);
            }

            // Barrier function, barr <- barr(x) where x o grad barr(x) = e
            static Real barr(Vector const & x) {
                return x.barr();
            }

            // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >=
            // 0} where y > 0.
            static Real srch(Vector const & x,Vector const & y) {
                return y.srch(x);
            }

            // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
            // operator.
            static void symm(Vector & x) {
                x.symm();
            }
        };

        // A couple of type shortcuts
        typedef Optizelle::ScalarValuedFunction <double,PythonVS>
            PyScalarValuedFunction;
        typedef Optizelle::VectorValuedFunction <double,PythonVS,PythonVS>
            PyVectorValuedFunction;
        typedef Optizelle::Operator <double,PythonVS,PythonVS> PyOperator;

        typedef Optizelle::Unconstrained <double,PythonVS>
            PyUnconstrained;
        typedef Optizelle::EqualityConstrained <double,PythonVS,PythonVS>
            PyEqualityConstrained;
        typedef Optizelle::InequalityConstrained <double,PythonVS,PythonVS>
            PyInequalityConstrained;
        typedef Optizelle::Constrained <double,PythonVS,PythonVS,PythonVS>
            PyConstrained;

        typedef Optizelle::json::Unconstrained <double,PythonVS>
            PyJsonUnconstrained;
        typedef Optizelle::json::EqualityConstrained <double,PythonVS,PythonVS>
            PyJsonEqualityConstrained;
        typedef Optizelle::json::InequalityConstrained<double,PythonVS,PythonVS>
            PyJsonInequalityConstrained;
        typedef Optizelle::json::Constrained <double,PythonVS,PythonVS,PythonVS>
            PyJsonConstrained;

        typedef typename Optizelle::RestartPackage <double>::t Reals;
        typedef typename Optizelle::RestartPackage <Natural>::t Naturals;
        typedef typename Optizelle::RestartPackage <std::string>::t Params;
        typedef typename Optizelle::RestartPackage <Vector>::t Vectors;

        // A simple scalar valued function interface, f : X -> R
        struct ScalarValuedFunction :
            public Optizelle::ScalarValuedFunction <double,PythonVS>
        {
        private:
            // Create some type shortcuts
            typedef Python::Vector Vector;

            // Data
            PyObjectPtr data;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(ScalarValuedFunction)

            // Create a function
            ScalarValuedFunction(PyObjectPtr const & data_);

            // <- f(x)
            double eval(Vector const & x) const;

            // g = grad f(x)
            void grad(Vector const & x,Vector & g) const;

            // H_dx = hess f(x) dx
            void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)const;
        };

        // A simple vector valued function interface, f : X -> Y
        struct VectorValuedFunction :
            public Optizelle::VectorValuedFunction<double,PythonVS,PythonVS>
        {
        private:
            // Create some type shortcuts
            typedef Python::Vector X_Vector;
            typedef Python::Vector Y_Vector;

            // Name of this function
            std::string const name;

            // Data
            PyObjectPtr data;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(VectorValuedFunction)

            // Create a function
            VectorValuedFunction(
                std::string const & name_,
                PyObjectPtr const & data_);

            // y=f(x)
            void eval(X_Vector const & x,Y_Vector & y) const;

            // y=f'(x)dx
            void p(X_Vector const & x,X_Vector const & dx,Y_Vector & y) const;

            // xhat=f'(x)*dy
            void ps(X_Vector const & x,const Y_Vector & dy,X_Vector & xhat)
                const;

            // xhat=(f''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                const Y_Vector & dy,
                X_Vector & xhat
            ) const;
        };

        // A linear operator specification, A : X->Y
        template <typename ProblemClass>
        struct Operator :
            public Optizelle::Operator <double,PythonVS,PythonVS>
        {
        private:
            // Create some type shortcuts
            typedef Python::Vector X_Vector;
            typedef Python::Vector Y_Vector;

            // Name of this function
            std::string const name;

            // Data
            PyObjectPtr data;

            // Optimization state.  Here's a funny trick.  Frequently, we
            // an operator like inv(g'(x)g'(x)*).  Notice, that this operator
            // depends on x, but the application of this operator doesn't
            // provide this information (Generally, the X_Vector is something
            // like dx.)  Hence, how do we get it?  In C++, we can just have a
            // reference to our variable hiding in the operator object.  When
            // interfacing to other languages, we can't.  In order to combat
            // this issue, we just pass the entire optimization state to the
            // operator application and then the user can extract what they
            // want.
            State <ProblemClass> & pystate;
            typename ProblemClass::State::t const & state;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Operator)

            // Create an operator
            Operator(
                std::string const & name_,
                PyObjectPtr const & data_,
                State <ProblemClass> & pystate_,
                typename ProblemClass::State::t const & state_
            ) :
                name(name_),
                data(data_),
                pystate(pystate_),
                state(state_)
            {}

            // y = A(x)
            void eval(X_Vector const & x,Y_Vector & y) const {
                // Convert the state to a Python state
                pystate.toPython(state);

                // Apply the operator to the state, x, and y
                auto eval = capi::PyObject_GetAttrString(data,"eval");
                capi::PyObject_CallObject3(
                    eval,
                    pystate.data,
                    x.data,
                    y.data,
                    __LOC__
                        + ", evaluation of the eval function in the operator "
                        + name + " failed");
            }
        };

        // Converts elements from C++ to Python
        namespace toPython {
            // Sets a real in a Python state
            void Real(
                std::string const & name,
                double const & value,
                PyObjectPtr & state
            );

            // Sets a natural in a Python state
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                PyObjectPtr & state
            );

            // Sets a parameter in a Python state
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<PyObjectPtr(enum_t const &)> const & toPython,
                enum_t const & value,
                PyObjectPtr & state
            ) {
                auto item = toPython(value);
                capi::PyObject_SetAttrString(state,name.c_str(),item);
            }

            // Sets a vector in a Python state
            void Vector(
                std::string const & name,
                Python::Vector const & value,
                PyObjectPtr & state
            );

            // Sets a list of vectors in a Python state
            void VectorList(
                std::string const & name,
                std::list <Python::Vector> const & values,
                PyObjectPtr & state
            );

            // Sets restart vectors in Python
            void Vectors(
                Python::Vectors const & values,
                PyObjectPtr & pyvalues
            );

            // Sets restart reals in Python
            void Reals(
                Python::Reals const & values,
                PyObjectPtr & pyvalues
            );

            // Sets restart naturals in Python
            void Naturals(
                Python::Naturals const & values,
                PyObjectPtr & pyvalues
            );

            // Sets restart parameters in Python
            void Params(
                Python::Params const & values,
                PyObjectPtr & pyvalues
            );
        }

        // Converts elements from Python to C++
        namespace fromPython {
            // Sets a real in a C++ state
            void Real(
                std::string const & name,
                PyObjectPtr const & pystate,
                double & value
            );

            // Sets a natural in a C++ state
            void Natural(
                std::string const & name,
                PyObjectPtr const & pystate,
                Optizelle::Natural & value
            );

            // Sets a param C++ state
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<enum_t(PyObjectPtr const &)> const & fromPython,
                PyObjectPtr const & pystate,
                enum_t & value
            ) {
                auto item = capi::PyObject_GetAttrString(pystate,name.c_str());
                value = fromPython(item);
            }

            // Sets a vector in a C++ state
            void Vector(
                std::string const & name,
                PyObjectPtr const & pystate,
                Python::Vector & value
            );

            // Sets a list of vectors in a C++ state
            void VectorList(
                std::string const & name,
                PyObjectPtr const & pystate,
                Python::Vector const & vec,
                std::list <Python::Vector> & values
            );

            // Sets a scalar-valued function in a C++ function bundle
            void ScalarValuedFunction(
                std::string const & name,
                PyObjectPtr const & pyfns,
                std::unique_ptr <PyScalarValuedFunction> & value
            );

            // Sets a vector-valued function in a C++ function bundle
            void VectorValuedFunction(
                std::string const & name,
                PyObjectPtr const & pyfns,
                std::unique_ptr <PyVectorValuedFunction> & value
            );

            // Sets an operator in a C++ function bundle
            template <typename ProblemClass>
            void Operator(
                std::string const & name,
                PyObjectPtr const & pyfns,
                State <ProblemClass> & pystate,
                typename ProblemClass::State::t const & state,
                std::unique_ptr <PyOperator> & value
            ) {
                value.reset(
                    new Python::Operator <ProblemClass> (
                        name,
                        capi::PyObject_GetAttrString(pyfns,name.c_str()),
                        pystate,
                        state));
            }

            // Sets restart vectors in C++
            void Vectors(
                Python::Vector const & vec,
                PyObjectPtr const & pyvalues,
                Python::Vectors & values
            );

            // Sets restart reals in C++
            void Reals(
                PyObjectPtr const & pyvalues,
                Python::Reals & values
            );

            // Sets restart naturals in C++
            void Naturals(
                PyObjectPtr const & pyvalues,
                Python::Naturals & values
            );

            // Sets restart parameters in C++
            void Params(
                PyObjectPtr const & pyvalues,
                Python::Params & values
            );
        }

        // Routines that manipulate and support problems of the form
        //
        // min_{x \in X} f(x)
        //
        // where f : X -> R.
        namespace Unconstrained {

            // Routines that manipulate the internal state of the optimization
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state
                void toPython_(
                    typename PyUnconstrained::State::t const & state,
                    PyObjectPtr & pystate);
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    Python::State <PyUnconstrained> & pystate);

                // Convert a Python state to C++
                void fromPython_(
                    PyObjectPtr const & pystate,
                    typename PyUnconstrained::State::t & state);
                void fromPython(
                    Python::State <PyUnconstrained> const & pystate,
                    typename PyUnconstrained::State::t & state);

                // Creates a state and inserts the elements into pystate
                PyObject * create(
                    PyObject * self,
                    PyObject * args);

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args);
            }

            // All the functions required by an optimization algorithm.
            namespace Functions{
                // Convert a Python bundle to C++
                template <typename ProblemClass>
                void fromPython_(
                    Python::Functions <ProblemClass> const & pyfns,
                    Python::State <ProblemClass> & pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyUnconstrained::Functions::t & fns
                ) {
                    fromPython::ScalarValuedFunction("f",pyfns.data,fns.f);
                    fromPython::Operator <ProblemClass> (
                        "PH",pyfns.data,pystate,state,fns.PH);
                }
                void fromPython(
                    Python::Functions <PyUnconstrained> const & pyfns,
                    Python::State <PyUnconstrained> & pystate,
                    typename PyUnconstrained::State::t const & state,
                    typename PyUnconstrained::Functions::t & fns
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                );
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                );

                // Capture data from structures controlled by the user.
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                );

                // Writes a json restart file
                PyObject* write_restart(
                    PyObject * self,
                    PyObject * args
                );

                // Reads a json restart file
                PyObject* read_restart(
                    PyObject * self,
                    PyObject * args
                );
            }
        }

        // Routines that manipulate and support problems of the form
        //
        // min_{x \in X} f(x) st g(x) = 0
        //
        // where f : X -> R and g : X -> Y
        namespace EqualityConstrained {
            // Routines that manipulate the internal state of the optimization
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state
                void toPython_(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObjectPtr & pystate);
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    Python::State <PyEqualityConstrained> & pystate);

                // Convert a Python state to C++
                void fromPython_(
                    PyObjectPtr const & pystate,
                    typename PyEqualityConstrained::State::t & state);
                void fromPython(
                    PyObjectPtr const & pystate,
                    typename PyEqualityConstrained::State::t & state);

                // Creates a state and inserts the elements into pystate
                PyObject * create(
                    PyObject * self,
                    PyObject * args);

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args);
            }

            // All the functions required by an optimization algorithm.
            namespace Functions{
                // Convert a Python bundle to C++
                template <typename ProblemClass>
                void fromPython_(
                    Python::Functions <ProblemClass> const & pyfns,
                    Python::State <ProblemClass> & pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyEqualityConstrained::Functions::t & fns
                ) {
                    fromPython::VectorValuedFunction("g",pyfns.data,fns.g);
                    fromPython::Operator <ProblemClass> ("PSchur_left",
                        pyfns.data,pystate,state,fns.PSchur_left);
                    fromPython::Operator <ProblemClass> ("PSchur_right",
                        pyfns.data,pystate,state,fns.PSchur_right);
                }
                void fromPython(
                    Python::Functions <PyEqualityConstrained> const & pyfns,
                    Python::State <PyEqualityConstrained> & pystate,
                    typename PyEqualityConstrained::State::t const & state,
                    typename PyEqualityConstrained::Functions::t & fns
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                );
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                );

                // Capture data from structures controlled by the user.
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                );

                // Writes a json restart file
                PyObject* write_restart(
                    PyObject * self,
                    PyObject * args
                );

                // Reads a json restart file
                PyObject* read_restart(
                    PyObject * self,
                    PyObject * args
                );
            }
        }

        // Routines that manipulate and support problems of the form
        //
        // min_{x \in X} f(x) st h(x) >=_K 0
        //
        // where f : X -> R and h : X -> Z
        namespace InequalityConstrained {
            // Routines that manipulate the internal state of the optimization
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state
                void toPython_(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObjectPtr & pystate);
                void toPython(
                    typename PyInequalityConstrained::State::t const & state,
                    Python::State <PyInequalityConstrained> & pystate);

                // Convert a Python state to C++
                void fromPython_(
                    PyObjectPtr const & pystate,
                    typename PyInequalityConstrained::State::t & state);
                void fromPython(
                    PyObjectPtr const & pystate,
                    typename PyInequalityConstrained::State::t & state);

                // Creates a state and inserts the elements into pystate
                PyObject * create(
                    PyObject * self,
                    PyObject * args);

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args);
            }

            // All the functions required by an optimization algorithm.
            namespace Functions{
                // Convert a Python bundle to C++
                template <typename ProblemClass>
                void fromPython_(
                    Python::Functions <ProblemClass> const & pyfns,
                    Python::State <ProblemClass> & pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyInequalityConstrained::Functions::t & fns
                ) {
                    fromPython::VectorValuedFunction("h",pyfns.data,fns.h);
                }
                void fromPython(
                    Python::Functions <PyInequalityConstrained> const & pyfns,
                    Python::State <PyInequalityConstrained> & pystate,
                    typename PyInequalityConstrained::State::t const & state,
                    typename PyInequalityConstrained::Functions::t & fns
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                );
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                );

                // Capture data from structures controlled by the user.
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                );

                // Writes a json restart file
                PyObject* write_restart(
                    PyObject * self,
                    PyObject * args
                );

                // Reads a json restart file
                PyObject* read_restart(
                    PyObject * self,
                    PyObject * args
                );
            }
        }

        // Routines that manipulate and support problems of the form
        // problem of the form
        //
        // min_{x \in X} f(x) st g(x) = 0, h(x) >=_K 0
        //
        // where f : X -> R, g : X -> Y, and h : X -> Z
        namespace Constrained {
            // Routines that manipulate the internal state of the optimization
            // algorithm
            namespace State {
                // Convert a C++ state to a Python state
                void toPython(
                    typename PyConstrained::State::t const & state,
                    Python::State <PyConstrained> & pystate);

                // Convert a Python state to C++
                void fromPython(
                    PyObjectPtr const & pystate,
                    typename PyConstrained::State::t & state);

                // Creates a state and inserts the elements into pystate
                PyObject * create(
                    PyObject * self,
                    PyObject * args);

                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args);
            }

            // All the functions required by an optimization algorithm.
            namespace Functions{
                // Convert a Python bundle to C++
                void fromPython(
                    Python::Functions <PyConstrained> const & pyfns,
                    Python::State <PyConstrained> & pystate,
                    typename PyConstrained::State::t const & state,
                    typename PyConstrained::Functions::t & fns
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                PyObject * getMin(
                    PyObject * self,
                    PyObject * args
                );
            }

            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user
                PyObject * release(
                    PyObject * self,
                    PyObject * args
                );

                // Capture data from structures controlled by the user.
                PyObject * capture(
                    PyObject * self,
                    PyObject * args
                );

                // Writes a json restart file
                PyObject* write_restart(
                    PyObject * self,
                    PyObject * args
                );

                // Reads a json restart file
                PyObject* read_restart(
                    PyObject * self,
                    PyObject * args
                );
            }
        }
    }
}
