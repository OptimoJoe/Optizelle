/*
Copyright 2013 OptimoJoe.

For the full copyright notice, see LICENSE.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Joseph Young (joe@optimojoe.com)
*/

#ifndef OPTIZELLE_PYTHON_H
#define OPTIZELLE_PYTHON_H

#include <Python.h>
#include <exception>
#include <optizelle/optizelle.h>
#include <optizelle/json.h>

// Alright, integrating C++ with Python is fraught with issues, but one that
// affects us specifically is const correctness.  There's not really a good
// way to handle this since Python doesn't have a concept of constant elements.
// In our code, we will consider something constant if we do not change its
// underlying value.  Now, we may attach to it, which increments its reference
// count, but we will consider this to be a constant operation.

namespace Optizelle {
    namespace StoppingCondition { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & opt_stop);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }

    namespace KrylovStop {
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & krylov_stop);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace KrylovSolverTruncated {
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & truncated_krylov);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }

    namespace AlgorithmClass { 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & algorithm_class);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace Operators{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & op);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace LineSearchDirection{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & dir);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace LineSearchKind{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & kind);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace OptimizationLocation{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & loc);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace InteriorPointMethod{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & ipm);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace CentralityStrategy{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & cstrat);

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace FunctionDiagnostics{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & diag); 

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }
    
    namespace DiagnosticScheme{ 
        // Converts t to a Python enumerated type
        PyObject * toPython(t const & dscheme); 

        // Converts a Python enumerated type to t 
        t fromPython(PyObject * const member);
    }

    namespace Python {
        
        // A function to alter the behavior of PyTuple_SetItem so that we don't
        // have to hand increment the reference to the object since SetItem
        // takes control of its arguments.
        void MyPyTuple_SetItem(PyObject * p,Natural const & pos,PyObject * o);

        // Like PyObject_GetAttrString, but returns obj.name1.name2 
        PyObject * PyObject_GetAttrString2(
            PyObject * const obj,
            std::string const & name1,
            std::string const & name2
        );

        // Calls a Python function with one argument 
        PyObject * PyObject_CallObject1(
            PyObject * const fn,
            PyObject * const arg1
        );
        
        // Calls a Python function with two arguments
        PyObject * PyObject_CallObject2(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2
        );
        
        // Calls a Python function with three arguments
        PyObject * PyObject_CallObject3(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2,
            PyObject * const arg3
        );
        
        // Calls a Python function with four arguments
        PyObject * PyObject_CallObject4(
            PyObject * const fn,
            PyObject * const arg1,
            PyObject * const arg2,
            PyObject * const arg3,
            PyObject * const arg4
        );

        // Used to catch Python exceptions
        struct Exception : public std::exception {
            Exception();
        };
        
        // Deep copy of a Python object and return the result
        PyObject * deepcopy(PyObject * const in);

        namespace PyObjectPtrMode {
            enum t : Natural {
                Capture,        // Capture the pointer
                Attach          // Attach to the pointer
            };
        }

        // A custom PyObject pointer that does proper clean-up on termination
        struct PyObjectPtr { 
        protected:
            // Internal storage of the pointer
            PyObject * ptr;
            
        public:
            // Disallow constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(PyObjectPtr);

            // On construction, initialize the pointer and figure out if
            // we're capturing the pointer or attaching to it
            explicit PyObjectPtr(
                PyObject * const ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            );

            // Move constructor
            explicit PyObjectPtr(PyObjectPtr&& ptr_) noexcept;

            // Move assignment operator
            PyObjectPtr const & operator = (PyObjectPtr&& ptr_) noexcept;

            // For a reset, we decrement the pointer and then assign a new
            // value.
            void reset(PyObject * const ptr_); 

            // For an attach, we decrement the pointer, assign a new value,
            // and then increment the reference count.
            void attach(PyObject * const ptr_);

            // On a get, we simply return the pointer.
            PyObject * get();
            
            // On a release, we return the underlying pointer and then clear
            // the vector.  This will prevent a decrement later.
            PyObject * release();

            // On destruction, decrement the Python reference counter and do
            // not delete the pointer.
            ~PyObjectPtr();
        };

        // A messaging utility that hooks directly into Python 
        struct Messaging : public Optizelle::Messaging, public PyObjectPtr {
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Messaging);

            // On construction, we just grab the pointer to the messaging object
            explicit Messaging(
                PyObject * const ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            );

            // Move constructor
            explicit Messaging(Messaging && msg) noexcept;

            // Move assignment operator
            Messaging const & operator = (Messaging && msg) noexcept;

            // Prints a message
            void print(std::string const & msg_) const;

            // Prints an error
            void error(std::string const & msg_) const;
        };

        // This class merges the vector space with a vector into a singular 
        // object.  We require this structure since Optizelle requires the
        // vector space to be static.  Since the user is passing us a vector
        // space dynamically, we merge the vector space functions with the
        // vectors and then statically define the vector space to call these
        // functions.
        struct Vector : public PyObjectPtr {
        private:
            // Messaging object
            Messaging msg;

            // Vector space
            PyObjectPtr vs;

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(Vector);

            // Create a vector with the appropriate messaging and vector space 
            explicit Vector(
                PyObject * const msg_,
                PyObject * const vs_,
                PyObject * const vec,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // Create a move constructor so we can interact with stl objects
            Vector(Vector && vec) noexcept;
            
            // Move assignment operator
            Vector const & operator = (Vector && vec) noexcept;

            // Memory allocation and size setting 
            Vector init();
            
            // y <- x (Shallow.  No memory allocation.)  Internal is y.
            void copy(Vector & x);

            // x <- alpha * x.  Internal is x.
            void scal(double const & alpha);

            // x <- 0.  Internal is x. 
            void zero();

            // y <- alpha * x + y.  Internal is y.
            void axpy(double const & alpha, Vector & x);

            // innr <- <x,y>.  Internal is y.
            double innr(Vector & x);

            // x <- random.  Internal is x.
            void rand();

            // Jordan product, z <- x o y.  Internal is z.
            void prod(Vector & x, Vector & y);

            // Identity element, x <- e such that x o e = x.  Internal is x.
            void id();

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y .
            // Internal is z.
            void linv(Vector & x, Vector & y); 

            // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
            // Internal is x.
            double barr();

            // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >=
            // 0} where y > 0.  Internal is y.
            double srch(Vector & x); 

            // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
            // operator.  Internal is x.
            void symm(); 
            
            // Converts (copies) a value into Python.  This assumes memory
            // has been allocated both in the vector as well as Python.
            void toPython(PyObject * const ptr);
            
            // Converts (copies) a value from Python.  This assumes memory
            // has been allocated both in the vector as well as Python.
            void fromPython(PyObject * const ptr);
        };
        
        // Python state
        template <typename ProblemClass>
        struct State : public PyObjectPtr {
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(State);

            // On construction, we just grab the pointer to the state object
            explicit State(
                PyObject * const ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            ) : PyObjectPtr(ptr_,mode) {}

            // Convert a C++ state to a Python state 
            void toPython(typename ProblemClass::State::t const & state);

            // Convert a Python state to C++ 
            void fromPython(typename ProblemClass::State::t & state);
        };
        
        // Python bundle of functions 
        template <typename ProblemClass>
        struct Functions : public PyObjectPtr {
        private:
            // Messaging object
            Messaging msg; 
            
            // Keep some states lying around so that we can communicate this
            // to our operator.
            State <ProblemClass> pystate;
            typename ProblemClass::State::t const & state;
            
        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Functions);

            // On construction, we just grab the pointer to the bundle object
            explicit Functions(
                PyObject * const msg_,
                PyObject * const pystate_,
                typename ProblemClass::State::t const & state_,
                PyObject * const ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            ) :
                PyObjectPtr(ptr_,mode),
                msg(msg_,PyObjectPtrMode::Attach),
                pystate(pystate_,PyObjectPtrMode::Attach),
                state(state_)
            {}

            // Convert a Python bundle to C++ 
            void fromPython(typename ProblemClass::Functions::t & fns);
        };

        // The state manipulator for Python
        template <typename ProblemClass>
        struct StateManipulator :
            public Optizelle::StateManipulator <ProblemClass>,
            public PyObjectPtr 
        {
        private:
            // Messaging object for reporting errors
            Messaging msg;

            // Keep a copy of a Python state lying around so that we can
            // use it to pass information back and forth to the Python
            // StateManipulator.
            mutable State <ProblemClass> pystate;

            // Similarly, we keep only the Python versin of the bundle of
            // functions lying around
            mutable PyObjectPtr pyfns;

        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(StateManipulator);

            // We need the Python state manipulator, a copy of a Python state
            // to pass information, and a copy of the Python functions.
            StateManipulator(
                PyObject * const msg_,
                PyObject * const pystate_,
                PyObject * const pyfns_,
                PyObject * const smanip_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            ) :
                PyObjectPtr(smanip_,mode),
                msg(msg_,PyObjectPtrMode::Attach),
                pystate(pystate_,PyObjectPtrMode::Attach),
                pyfns(pyfns_,PyObjectPtrMode::Attach)
            {}

            // Application
            void eval(
                const typename ProblemClass::Functions::t & fns,
                typename ProblemClass::State::t & state,
                OptimizationLocation::t const & loc_
            ) const {
                // Convert the C++ state to a Python state
                pystate.toPython(state);

                // Convert the lcoation to Python
                PyObjectPtr loc(OptimizationLocation::toPython(loc_));
            
                // Call the Python state manipulator 
                // give it pystate and pyfns.  Note, pyfns is given raw.
                PyObjectPtr eval(PyObject_GetAttrString(ptr,"eval"));
                PyObjectPtr ret(PyObject_CallObject3(
                    eval.get(),
                    pyfns.get(),
                    pystate.get(),
                    loc.get()));

                // Check errors
                if(ret.get()==nullptr)
                    msg.error("Evaluation of the StateManipulator object "
                        "failed.");

                // Convert the Python state to the C++ state 
                pystate.fromPython(state);
            }
        };

        // Vector space that works with Python objects
        template <typename Real=double> 
        struct PythonVS { 
            // Prevent constructors 
            NO_CONSTRUCTORS(PythonVS);

            // Setup the vector 
            typedef Optizelle::Python::Vector Vector; 

            // Memory allocation and size setting 
            static Vector init(Vector const & x) { 
                return std::move(const_cast <Vector &> (x).init());
            } 

            // y <- x (Shallow.  No memory allocation.) 
            static void copy(Vector const & x, Vector & y) { 
                y.copy(const_cast <Vector &> (x));
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
                y.axpy(alpha,const_cast <Vector &> (x));
            } 

            // innr <- <x,y> 
            static Real innr(Vector const & x,Vector const & y) { 
                return const_cast <Vector &> (y).innr(const_cast <Vector &>(x));
            } 
        
            // x <- random
            static void rand(Vector & x){
                x.rand();
            }

            // Jordan product, z <- x o y 
            static void prod(Vector const & x, Vector const & y, Vector & z) { 
                z.prod(const_cast <Vector &> (x),const_cast <Vector &> (y));
            } 

            // Identity element, x <- e such that x o e = x 
            static void id(Vector & x) { 
                x.id();
            } 

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
            static void linv(Vector const & x, Vector const & y, Vector & z) { 
                z.linv(const_cast <Vector &> (x),const_cast <Vector &> (y));
            } 

            // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
            static Real barr(Vector const & x) { 
                return const_cast <Vector &> (x).barr();
            } 

            // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >=
            // 0} where y > 0. 
            static Real srch(Vector const & x,Vector const & y) {  
                return const_cast <Vector &> (y).srch(const_cast <Vector &>(x));
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
            public Optizelle::ScalarValuedFunction <double,PythonVS>,
            public PyObjectPtr
        {
        private:
            // Create some type shortcuts
            typedef PythonVS <>::Vector Vector; 

            // Messaging object
            Messaging msg; 

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(ScalarValuedFunction);

            // Create a function 
            explicit ScalarValuedFunction(
                PyObject * const msg_,
                PyObject * const f,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // <- f(x) 
            double eval(Vector const & x) const; 

            // g = grad f(x) 
            void grad(Vector const & x,Vector & g) const; 

            // H_dx = hess f(x) dx 
            void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)const;
        };

        // A simple vector valued function interface, f : X -> Y
        struct VectorValuedFunction :
            public Optizelle::VectorValuedFunction<double,PythonVS,PythonVS>,
            public PyObjectPtr
        {
        private:
            // Create some type shortcuts
            typedef PythonVS <>::Vector X_Vector; 
            typedef PythonVS <>::Vector Y_Vector; 

            // Messaging object
            Messaging msg; 

            // Name of this function
            std::string const name;

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(VectorValuedFunction);

            // Create a function 
            explicit VectorValuedFunction(
                std::string const & name_,
                PyObject * const msg_,
                PyObject * const f,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // y=f(x)
            void eval(X_Vector const & x,Y_Vector & y) const;

            // y=f'(x)dx 
            void p(X_Vector const & x,X_Vector const & dx,Y_Vector & y) const;

            // z=f'(x)*dy
            void ps(X_Vector const & x,const Y_Vector & dy,X_Vector & z) const; 
             
            // z=(f''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                const Y_Vector & dy,
                X_Vector & z
            ) const; 
        };
        
        // A linear operator specification, A : X->Y 
        template <typename ProblemClass>
        struct Operator :
            public Optizelle::Operator <double,PythonVS,PythonVS>,
            public PyObjectPtr
        {
        private:
            // Create some type shortcuts
            typedef PythonVS <>::Vector X_Vector; 
            typedef PythonVS <>::Vector Y_Vector; 

            // Messaging object
            Messaging msg; 

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
            mutable State <ProblemClass> pystate;
            typename ProblemClass::State::t const & state;

            // Name of this function
            std::string const name;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Operator);

            // Create an operator 
            explicit Operator(
                std::string const & name_,
                PyObject * const msg_,
                PyObject * const op,
                PyObject * const pystate_,
                typename ProblemClass::State::t const & state_,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture
            ) :
                PyObjectPtr(op,mode),
                msg(msg_,PyObjectPtrMode::Attach),
                pystate(pystate_,PyObjectPtrMode::Attach),
                state(state_),
                name(name_)
            {}

            // y = A(x)
            void eval(X_Vector const & x,Y_Vector & y) const {
                // Convert the state to a Python state
                pystate.toPython(state);

                // Apply the operator to the state, x, and y
                PyObjectPtr eval(PyObject_GetAttrString(ptr,"eval"));
                PyObjectPtr ret(PyObject_CallObject3(
                    eval.get(),
                    pystate.get(),
                    const_cast <X_Vector &> (x).get(),
                    y.get()));

                // Check errors
                if(ret.get()==nullptr) {
                    std::stringstream ss;
                    ss << "Evaluation of the eval function in the operator "
                        << name << " failed.";
                    msg.error(ss.str());
                }
            }
        };
        
        // Calls the Optizelle exception with a string
        void PyErr_SetString_Optizelle(std::string const & msg);

        // Converts an Optizelle enumerated type to a PyObject * 
        PyObject * enumToPyObject(
            std::string const & type,
            std::string const & member 
        );
        
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        );

        // Converts elements from C++ to Python 
        namespace toPython {
            // Sets a real in a Python state 
            void Real(
                std::string const & name,
                double const & value,
                PyObject * const obj 
            );

            // Sets a natural in a Python state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                PyObject * const obj 
            );
           
            // Sets a parameter in a Python state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<PyObject *(enum_t const &)> const & toPython,
                enum_t const & value,
                PyObject * const obj
            ) {
                PyObjectPtr item(toPython(value));
                PyObject_SetAttrString(obj,name.c_str(),item.get());
            }
            
            // Sets a vector in a Python state 
            void Vector(
                std::string const & name,
                Python::Vector const & value,
                PyObject * const obj 
            );
            
            // Sets a list of vectors in a Python state 
            void VectorList(
                std::string const & name,
                std::list <Python::Vector> const & values,
                PyObject * const obj 
            );
            
            // Sets a scalar-valued function in a Python function bundle 
            void ScalarValuedFunction(
                std::string const & name,
                PyObject * const msg,
                PyObject * const obj,
                std::unique_ptr <PyScalarValuedFunction> & value
            );
            
            // Sets a vector-valued function in a Python function bundle 
            void VectorValuedFunction(
                std::string const & name,
                PyObject * const msg,
                PyObject * const obj,
                std::unique_ptr <PyVectorValuedFunction> & value
            );
            
            // Sets an operator in a Python function bundle 
            template <typename ProblemClass>
            void Operator(
                std::string const & name,
                PyObject * const msg,
                PyObject * const obj,
                PyObject * const pystate,
                typename ProblemClass::State::t const & state,
                std::unique_ptr <PyOperator> & value
            ) {
                value.reset(new Python::Operator <ProblemClass> (
                    name,
                    msg,
                    PyObject_GetAttrString(obj,name.c_str()),
                    pystate,
                    state));
            }
        
            // Sets restart vectors in Python 
            void Vectors(
                Python::Vectors const & values,
                PyObject * const pyvalues 
            );

            // Sets restart reals in Python 
            void Reals(
                Python::Reals const & values,
                PyObject * const pyvalues 
            );
            
            // Sets restart naturals in Python 
            void Naturals(
                Python::Naturals const & values,
                PyObject * const pyvalues 
            );
            
            // Sets restart parameters in Python 
            void Params(
                Python::Params const & values,
                PyObject * const pyvalues 
            );
        }

        // Converts elements from Python to C++ 
        namespace fromPython {
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                PyObject * const obj,
                double & value
            );
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                PyObject * const obj,
                Optizelle::Natural & value
            );
           
            // Sets a param C++ state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<enum_t(PyObject * const)> const & fromPython,
                PyObject * const obj,
                enum_t & value
            ) {
                PyObjectPtr item(PyObject_GetAttrString(
                    const_cast <PyObject *> (obj),name.c_str()));
                value = fromPython(item.get());
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                PyObject * const obj,
                Python::Vector & value
            );
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                PyObject * const obj,
                Python::Vector const & vec,
                std::list <Python::Vector> & values
            );
        
            // Sets restart vectors in C++ 
            void Vectors(
                Python::Vector const & vec,
                PyObject * const pyvalues,
                Python::Vectors & values
            );
            
            // Sets restart reals in C++ 
            void Reals(
                PyObject * const pyvalues,
                Python::Reals & values
            );
            
            // Sets restart naturals in C++ 
            void Naturals(
                PyObject * const pyvalues,
                Python::Naturals & values
            );
            
            // Sets restart parameters in C++ 
            void Params(
                PyObject * const pyvalues,
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
                    PyObject * const pystate
                );
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    PyObject * const pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyUnconstrained::State::t & state
                );
                void fromPython(
                    PyObject * const pystate,
                    typename PyUnconstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                );
                
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                template <typename ProblemClass>
                void fromPython_(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyUnconstrained::Functions::t & fns 
                ) {
                    toPython::ScalarValuedFunction("f",msg,pyfns,fns.f);
                    toPython::Operator <ProblemClass> (
                        "PH",msg,pyfns,pystate,state,fns.PH);
                }
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
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
                    PyObject * const pystate
                );
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObject * const pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyEqualityConstrained::State::t & state
                );
                void fromPython(
                    PyObject * const pystate,
                    typename PyEqualityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                );
                
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                template <typename ProblemClass>
                void fromPython_(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyEqualityConstrained::Functions::t & fns 
                ) {
                    toPython::VectorValuedFunction("g",msg,pyfns,fns.g);
                    toPython::Operator <ProblemClass> ("PSchur_left",
                        msg,pyfns,pystate,state,fns.PSchur_left);
                    toPython::Operator <ProblemClass> ("PSchur_right",
                        msg,pyfns,pystate,state,fns.PSchur_right);
                }
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
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
                    PyObject * const pystate
                );
                void toPython(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObject * const pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject * const pystate,
                    typename PyInequalityConstrained::State::t & state
                );
                void fromPython(
                    PyObject * const pystate,
                    typename PyInequalityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                );
                
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                template <typename ProblemClass>
                void fromPython_(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
                    typename ProblemClass::State::t const & state,
                    typename PyInequalityConstrained::Functions::t & fns 
                ) {
                    toPython::VectorValuedFunction("h",msg,pyfns,fns.h);
                }
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
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
                    PyObject * const pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython(
                    PyObject * const pystate,
                    typename PyConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject * create(
                    PyObject * self,
                    PyObject * args
                );
                
                // Read json parameters from file
                PyObject * readJson(
                    PyObject * self,
                    PyObject * args
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Python bundle to C++ 
                void fromPython(
                    PyObject * const msg,
                    PyObject * const pyfns,
                    PyObject * const pystate,
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
#endif
