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
        
// A macro to alter the behavior of PyTuple_SetItem so that we don't
// have to hand increment the reference to the object since SetItem
// takes control of its arguments.
#define MyPyTuple_SetItem(p,pos,o) \
        Py_INCREF(o); \
        PyTuple_SetItem(p,pos,o);

namespace Optizelle {
    namespace StoppingCondition { 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython : public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }

    namespace KrylovStop {
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython : public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace KrylovSolverTruncated {
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython : public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }

    namespace AlgorithmClass { 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace Operators{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace LineSearchDirection{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace LineSearchKind{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace OptimizationLocation{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & loc) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace InteriorPointMethod{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }
    
    namespace CentralityStrategy{ 
        // Converts t to a Python enumerated type
        struct toPython : public std::unary_function <t const &,PyObject *> {
            PyObject* operator () (t const & opt_stop) const;
        };

        // Converts a Python enumerated type to t 
        struct fromPython: public std::unary_function
            <PyObject const * const,t>
        {
            t operator () (PyObject const * const member) const;
        };
    }

    namespace Python {
        // Like PyObject_GetAttrString, but returns obj.name1.name2 
        PyObject* PyObject_GetAttrString2(
            PyObject const * const obj,
            std::string const & name1,
            std::string const & name2
        );

        // Calls a Python function with one argument 
        PyObject* PyObject_CallObject1(
            PyObject const * const fn,
            PyObject const * const arg1
        );
        
        // Calls a Python function with two arguments
        PyObject* PyObject_CallObject2(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2
        );
        
        // Calls a Python function with three arguments
        PyObject* PyObject_CallObject3(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2,
            PyObject const * const arg3
        );
        
        // Calls a Python function with four arguments
        PyObject* PyObject_CallObject4(
            PyObject const * const fn,
            PyObject const * const arg1,
            PyObject const * const arg2,
            PyObject const * const arg3,
            PyObject const * const arg4
        );

        // Used to catch Python exceptions
        struct Exception : public std::exception {
            Exception();
        };
        
        // Deep copy of a Python object and return the result
        PyObject* deepcopy(PyObject const * const in);

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
            PyObject* ptr;
            
        public:
            // Disallow constructors 
            PyObjectPtr() = delete;
            PyObjectPtr(PyObjectPtr const &) = delete;
            void operator = (PyObjectPtr const &) = delete;

            // On construction, initialize the pointer and figure out if
            // we're capturing the pointer or attaching to it
            explicit PyObjectPtr(
                PyObject* ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            );

            // Move constructor
            explicit PyObjectPtr(PyObjectPtr&& ptr_) noexcept;

            // For a reset, we decrement the pointer and then assign a new
            // value.
            void reset(PyObject* ptr_); 

            // For an attach, we decrement the pointer, assign a new value,
            // and then increment the reference count.
            void attach(PyObject* ptr_);

            // On a get, we simply return the pointer.
            PyObject const * const get() const;
            PyObject* get();
            
            // On a release, we return the underlying pointer and then clear
            // the vector.  This will prevent a decrement later.
            PyObject* release();

            // On destruction, decrement the Python reference counter and do
            // not delete the pointer.
            ~PyObjectPtr();
        };

        // A messaging utility that hooks directly into Python 
        struct Messaging : public Optizelle::Messaging, public PyObjectPtr {
            // Disallow constructors
            Messaging() = delete;
            Messaging(Messaging const &) = delete;
            void operator = (Messaging const &) = delete;

            // On construction, we just grab the pointer to the messaging object
            explicit Messaging(
                PyObject* ptr_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            );

            // Move constructor
            explicit Messaging(Messaging&& ptr_) noexcept;

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
            Vector(Vector const &) = delete;
            Vector& operator = (Vector const &) = delete;

            // Create an empty vector.  Don't use these objects until
            // initialized.
            explicit Vector();

            // Create a vector with the appropriate messaging and vector space 
            explicit Vector(PyObject* msg_,PyObject* vs_,PyObject* vec,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // Create a move constructor so we can interact with stl objects
            explicit Vector(Vector&& vec);

            // Memory allocation and size setting 
            void init(const Vector& vec);
            
            // y <- x (Shallow.  No memory allocation.) 
            void copy(const Vector& x);

            // x <- alpha * x
            void scal(const double& alpha);

            // x <- 0 
            void zero();

            // y <- alpha * x + y 
            void axpy(const double& alpha, const Vector& x);

            // innr <- <x,y> 
            double innr(const Vector& x) const;

            // Jordan product, z <- x o y
            void prod(const Vector& x, const Vector& y);

            // Identity element, x <- e such that x o e = x 
            void id();

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
            void linv(const Vector& x, const Vector& y); 

            // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
            double barr() const;

            // Line search, srch <- argmax {alpha in Real >= 0 : alpha x + y >=
            // 0} where y > 0.  
            double srch(const Vector& x) const; 

            // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
            // operator.
            void symm(); 
        };
        
        // Python state
        template <typename ProblemClass>
        struct State : public PyObjectPtr {
            // Disallow constructors
            State() = delete;
            State(State const &) = delete;
            void operator = (State const &) = delete;

            // On construction, we just grab the pointer to the state object
            explicit State(
                PyObject* ptr_,
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
            Functions() = delete;
            Functions(Functions const &) = delete;
            void operator = (Functions const &) = delete;

            // On construction, we just grab the pointer to the bundle object
            explicit Functions(
                PyObject* msg_,
                PyObject* pystate_,
                typename ProblemClass::State::t const & state_,
                PyObject* ptr_,
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
            PyObjectPtr pyfns;

        public:
            // Disallow constructors
            StateManipulator() = delete;
            StateManipulator(StateManipulator const &) = delete;
            void operator = (StateManipulator const &) = delete;

            // We need the Python state manipulator, a copy of a Python state
            // to pass information, and a copy of the Python functions.
            StateManipulator(
                PyObject* msg_,
                PyObject* pystate_,
                PyObject* pyfns_,
                PyObject* smanip_,
                PyObjectPtrMode::t const mode = PyObjectPtrMode::Capture
            ) :
                PyObjectPtr(smanip_,mode),
                msg(msg_,PyObjectPtrMode::Attach),
                pystate(pystate_,PyObjectPtrMode::Attach),
                pyfns(pyfns_,PyObjectPtrMode::Attach)
            {}

            // Application
            void operator () (
                const typename ProblemClass::Functions::t& fns,
                typename ProblemClass::State::t& state,
                OptimizationLocation::t loc
            ) const {
                // Convert the C++ state to a Python state
                pystate.toPython(state);

                // Convert the lcoation to Python
                PyObjectPtr loc_(OptimizationLocation::toPython()(loc));
            
                // Call the Python state manipulator 
                // give it pystate and pyfns.  Note, pyfns is given raw.
                PyObjectPtr eval_(PyObject_GetAttrString(ptr,"eval"));
                PyObjectPtr ret(PyObject_CallObject3(
                    eval_.get(),
                    pyfns.get(),
                    pystate.get(),
                    loc_.get()));

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
            PythonVS() = delete;
            PythonVS(PythonVS const &) = delete;
            void operator = (PythonVS const &) = delete;

            // Setup the vector 
            typedef Optizelle::Python::Vector Vector; 

            // Memory allocation and size setting 
            static void init(Vector const & x, Vector& y) { 
                y.init(x);
            } 

            // y <- x (Shallow.  No memory allocation.) 
            static void copy(Vector const & x, Vector& y) { 
                y.copy(x);
            } 

            // x <- alpha * x 
            static void scal(Real const & alpha, Vector& x) { 
                x.scal(alpha);
            } 

            // x <- 0 
            static void zero(Vector& x) { 
                x.zero();
            } 

            // y <- alpha * x + y 
            static void axpy(Real const & alpha, Vector const & x, Vector& y) { 
                y.axpy(alpha,x);
            } 

            // innr <- <x,y> 
            static Real innr(Vector const & x,Vector const & y) { 
                return y.innr(x);
            } 

            // Jordan product, z <- x o y 
            static void prod(Vector const & x, Vector const & y, Vector& z) { 
                z.prod(x,y);
            } 

            // Identity element, x <- e such that x o e = x 
            static void id(Vector& x) { 
                x.id();
            } 

            // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
            static void linv(Vector const & x, Vector const & y, Vector& z) { 
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
            static void symm(Vector& x) {
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
            ScalarValuedFunction() = delete;
            ScalarValuedFunction(ScalarValuedFunction const &) = delete;
            void operator = (ScalarValuedFunction const &) = delete;

            // Create a function 
            explicit ScalarValuedFunction(
                PyObject* msg_,
                PyObject* f,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // <- f(x) 
            double operator () (const Vector& x) const; 

            // g = grad f(x) 
            void grad(const Vector& x,Vector& g) const; 

            // H_dx = hess f(x) dx 
            void hessvec(const Vector& x,const Vector& dx,Vector& H_dx) const;
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
            VectorValuedFunction() = delete;
            VectorValuedFunction(VectorValuedFunction const &) = delete;
            void operator = (VectorValuedFunction const &) = delete;

            // Create a function 
            explicit VectorValuedFunction(
                std::string const & name_,
                PyObject* msg_,
                PyObject* f,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture);

            // y=f(x)
            void operator () (const X_Vector& x,Y_Vector& y) const ;

            // y=f'(x)dx 
            void p(const X_Vector& x,const X_Vector& dx,Y_Vector& y) const;

            // z=f'(x)*dy
            void ps(const X_Vector& x,const Y_Vector& dy,X_Vector& z)
                const; 
             
            // z=(f''(x)dx)*dy
            void pps(
                const X_Vector& x,
                const X_Vector& dx,
                const Y_Vector& dy,
                X_Vector& z
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
            Operator() = delete;
            Operator(Operator const &) = delete;
            void operator = (Operator const &) = delete;

            // Create an operator 
            explicit Operator(
                std::string const & name_,
                PyObject* msg_,
                PyObject* op,
                PyObject* pystate_,
                typename ProblemClass::State::t const & state_,
                PyObjectPtrMode::t mode=PyObjectPtrMode::Capture
            ) :
                PyObjectPtr(op,mode),
                msg(msg_,PyObjectPtrMode::Attach),
                pystate(pystate_,PyObjectPtrMode::Attach),
                state(state_),
                name(name_)
            { }

            // y = A(x)
            void operator () (const X_Vector& x,Y_Vector &y) const {
                // Convert the state to a Python state
                pystate.toPython(state);

                // Apply the operator to the state, x, and y
                PyObjectPtr eval_(PyObject_GetAttrString(ptr,"eval"));
                PyObjectPtr ret(PyObject_CallObject3(
                    eval_.get(),pystate.get(),x.get(),y.get()));

                // Check errors
                if(ret.get()==nullptr) {
                    std::stringstream ss;
                    ss << "Evaluation of the eval function of the operator "
                        << name << "failed.";
                    msg.error(ss.str());
                }
            }
        };
        
        // Calls the Optizelle exception with a string
        void PyErr_SetString_Optizelle(std::string const& msg);

        // Converts an Optizelle enumerated type to a PyObject* 
        PyObject* enumToPyObject(
            std::string const & type,
            std::string const & member 
        );
        
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        ); 

        // Sets a floating point in a Python class 
        void setFloat(
            std::string const & name,
            double const & value,
            PyObject* obj 
        );
        
        // Sets a floating point in a C++ class 
        void setFloat(
            std::string const & name,
            PyObject const * const obj,
            double & value
        );
        
        // Sets an integer in a Python class 
        void setNatural(
            std::string const & name,
            Natural const & value,
            PyObject* obj 
        );
        
        // Sets an integer in a C++ class 
        void setNatural(
            std::string const & name,
            PyObject const * const obj,
            Natural & value
        );
       
        // Sets an enumerated value in a Python class
        template <typename enum_t>
        void setEnum(
            std::string const & name,
            std::function<PyObject*(enum_t const &)> const & toPython,
            enum_t const & value,
            PyObject *obj
        ) {
            PyObjectPtr item(toPython(value));
            PyObject_SetAttrString(obj,name.c_str(),item.get());
        }
       
        // Sets an enumerated value in a C++ class
        template <typename enum_t>
        void setEnum(
            std::string const & name,
            std::function<enum_t(PyObject const * const)> const & fromPython,
            PyObject const * const obj,
            enum_t & value
        ) {
            PyObjectPtr item(PyObject_GetAttrString(
                const_cast <PyObject*> (obj),name.c_str()));
            value = fromPython(item.get());
        }
        
        // Sets a vector in a Python class 
        void setVector(
            std::string const & name,
            Vector const & value,
            PyObject* obj 
        );
        
        // Sets a vector in a C++ class 
        void setVector(
            std::string const & name,
            PyObject const * const obj,
            Vector & value
        );
        
        // Sets a list of vectors in a Python class 
        void setVectors(
            std::string const & name,
            std::list <Vector> const & values,
            PyObject* obj 
        );
        
        // Sets a list of vectors in a C++ class 
        void setVectors(
            std::string const & name,
            PyObject const * const obj,
            std::list <Vector> & values
        );
        
        // Sets a scalar-valued function in a Python class 
        void setScalarValuedFunction(
            std::string const & name,
            PyObject * const msg,
            PyObject * const obj,
            std::unique_ptr <PyScalarValuedFunction> & value
        );
        
        // Sets a vector-valued function in a Python class 
        void setVectorValuedFunction(
            std::string const & name,
            PyObject * const msg,
            PyObject * const obj,
            std::unique_ptr <PyVectorValuedFunction> & value
        );
        
        // Sets a linear operator in a Python class 
        template <typename ProblemClass>
        void setOperator(
            std::string const & name,
            PyObject * const msg,
            PyObject * const obj,
            PyObject * const pystate,
            typename ProblemClass::State::t const & state,
            std::unique_ptr <PyOperator> & value
        ) {
            value.reset(new Operator <ProblemClass> (
                name,
                msg,
                PyObject_GetAttrString(obj,name.c_str()),
                pystate,
                state));
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
                    PyObject* pystate
                );
                void toPython(
                    typename PyUnconstrained::State::t const & state,
                    PyObject* pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyUnconstrained::State::t & state
                );
                void fromPython(
                    PyObject const * const pystate,
                    typename PyUnconstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
                );
                
                // Read json parameters from file
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                    setScalarValuedFunction("f",msg,pyfns,fns.f);
                    setOperator <ProblemClass> (
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
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
                    PyObject* pystate
                );
                void toPython(
                    typename PyEqualityConstrained::State::t const & state,
                    PyObject* pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyEqualityConstrained::State::t & state
                );
                void fromPython(
                    PyObject const * const pystate,
                    typename PyEqualityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
                );
                
                // Read json parameters from file
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                    setVectorValuedFunction("g",msg,pyfns,fns.g);
                    setOperator <ProblemClass> ("PSchur_left",
                        msg,pyfns,pystate,state,fns.PSchur_left);
                    setOperator <ProblemClass> ("PSchur_right",
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
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
                    PyObject* pystate
                );
                void toPython(
                    typename PyInequalityConstrained::State::t const & state,
                    PyObject* pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython_(
                    PyObject const * const pystate,
                    typename PyInequalityConstrained::State::t & state
                );
                void fromPython(
                    PyObject const * const pystate,
                    typename PyInequalityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
                );
                
                // Read json parameters from file
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                    setVectorValuedFunction("h",msg,pyfns,fns.h);
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
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
                    PyObject* pystate
                );
                
                // Convert a Python state to C++ 
                void fromPython(
                    PyObject const * const pystate,
                    typename PyConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into pystate 
                PyObject* create(
                    PyObject* self,
                    PyObject* args
                );
                
                // Read json parameters from file
                PyObject* readJson(
                    PyObject* self,
                    PyObject* args
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
                PyObject* getMin(
                    PyObject* self,
                    PyObject* args
                );
            }
        }
    }
}
#endif
