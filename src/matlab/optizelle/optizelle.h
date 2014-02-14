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

#ifndef OPTIZELLE_MATLAB_H
#define OPTIZELLE_MATLAB_H

#include "mex.h"
#include "optizelle/optizelle.h"
#include "optizelle/json.h"

namespace Optizelle {
    namespace StoppingCondition { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & opt_stop);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }

    namespace KrylovStop {
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & krylov_stop);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace KrylovSolverTruncated {
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & truncated_krylov);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }

    namespace AlgorithmClass { 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & algorithm_class);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace Operators{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & op);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace LineSearchDirection{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & dir);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace LineSearchKind{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & kind);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace OptimizationLocation{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & loc);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace InteriorPointMethod{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & ipm);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace CentralityStrategy{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & cstrat);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace FunctionDiagnostics{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & diag); 

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }
    
    namespace DiagnosticScheme{ 
        // Converts t to a Matlab enumerated type
        mxArray * toMatlab(t const & dscheme); 

        // Converts a Matlab enumerated type to t 
        t fromMatlab(mxArray * const member);
    }

    namespace Matlab {
        // Used to catch Matlab exceptions
        struct Exception : public std::exception {
            Exception();
        };

        // Calls a Matlab function with one argument 
        std::pair <mxArray *,int> mxArray_CallObject1(
            mxArray * const fn,
            mxArray * const arg1
        );
        
        // Calls a Matlab function with two arguments
        std::pair <mxArray *,int> mxArray_CallObject2(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2
        );
        
        // Calls a Matlab function with three arguments
        std::pair <mxArray *,int> mxArray_CallObject3(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3
        );
        
        // Calls a Matlab function with four arguments
        std::pair <mxArray *,int> mxArray_CallObject4(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3,
            mxArray * const arg4
        );

        // Creates a Matlab double from a C++ double
        mxArray * mxArray_FromDouble(double const x_);

        // Creates a Matlab int from a C++ size_t 
        mxArray * mxArray_FromSize_t(Natural const x_);

        // Imports a piece of the Optizelle module.  This makes a deep
        // copy of the eventual imported object, so the result needs
        // to have its memory managed. 
        mxArray * importOptizelle(std::string const & module);

        // Converts an Optizelle enumerated type to a mxArray * 
        mxArray * enumToMxArray(
            std::string const & type,
            std::string const & member 
        );
        
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        );

        namespace mxArrayPtrMode {
            enum t : Natural {
                Capture,        // Capture the pointer
                Attach          // Attach to the pointer
            };
        }

        // A custom mxArray pointer that does proper clean-up on termination
        struct mxArrayPtr { 
        protected:
            // Internal storage of the pointer
            mxArray * ptr;

            // Whether or not we're attached to or controlling this pointer
            mxArrayPtrMode::t mode;
            
        public:
            // Disallow constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(mxArrayPtr);

            // On construction, initialize the pointer and figure out if
            // we're capturing the pointer or attaching to it
            mxArrayPtr(
                mxArray * const ptr_,
                mxArrayPtrMode::t const mode = mxArrayPtrMode::Capture
            );

            // Move constructor
            mxArrayPtr(mxArrayPtr&& ptr_) noexcept;

            // Move assignment operator
            mxArrayPtr const & operator = (mxArrayPtr&& ptr_) noexcept;

            // For a reset, we destroy the pointer and then assign a new
            // value.
            void reset(mxArray * const ptr_); 

            // For an attach, we destroy the pointer then assign a new value
            void attach(mxArray * const ptr_);

            // On a get, we simply return the pointer.
            mxArray * get();
            
            // On a release, we return the underlying pointer and then clear
            // the vector.  This will prevent destruction later. 
            mxArray * release();

            // On destruction, destroy the pointer. 
            ~mxArrayPtr();
        };

        // A messaging utility that hooks directly into Matlab 
        struct Messaging : public Optizelle::Messaging, public mxArrayPtr {
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Messaging);

            // On construction, we just grab the pointer to the messaging object
            explicit Messaging(
                mxArray * const ptr_,
                mxArrayPtrMode::t const mode = mxArrayPtrMode::Capture
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
        struct Vector : public mxArrayPtr {
        private:
            // Messaging object
            Messaging msg;

            // Vector space
            mxArrayPtr vs;

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(Vector);

            // Create a vector with the appropriate messaging and vector space 
            explicit Vector(
                mxArray * const msg_,
                mxArray * const vs_,
                mxArray * const vec,
                mxArrayPtrMode::t mode=mxArrayPtrMode::Capture);

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
            
            // Converts (copies) a value into Matlab.  
            mxArray * toMatlab();
            
            // Converts (copies) a value from Matlab.  This assumes that the
            // vector space functions have already been properly assigned.
            void fromMatlab(mxArray * const ptr);
        };
        
        // Matlab state
        template <typename ProblemClass>
        struct State : public mxArrayPtr {
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(State);

            // On construction, we just grab the pointer to the state object
            explicit State(
                mxArray * const ptr_,
                mxArrayPtrMode::t const mode = mxArrayPtrMode::Capture
            ) : mxArrayPtr(ptr_,mode) {}

            // Convert a C++ state to a Matlab state 
            void toMatlab(typename ProblemClass::State::t const & state);

            // Convert a Matlab state to C++ 
            void fromMatlab(typename ProblemClass::State::t & state);
        };
        
        // Matlab bundle of functions 
        template <typename ProblemClass>
        struct Functions : public mxArrayPtr {
        private:
            // Messaging object
            Messaging msg; 
            
            // Keep some states lying around so that we can communicate this
            // to our operator.
            State <ProblemClass> mxstate;
            typename ProblemClass::State::t const & state;
            
        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Functions);

            // On construction, we just grab the pointer to the bundle object
            explicit Functions(
                mxArray * const msg_,
                mxArray * const mxstate_,
                typename ProblemClass::State::t const & state_,
                mxArray * const ptr_,
                mxArrayPtrMode::t const mode = mxArrayPtrMode::Capture
            ) :
                mxArrayPtr(ptr_,mode),
                msg(msg_,mxArrayPtrMode::Attach),
                mxstate(mxstate_,mxArrayPtrMode::Attach),
                state(state_)
            {}

            // Convert a Matlab bundle to C++ 
            void fromMatlab(typename ProblemClass::Functions::t & fns);
        };
        
        // The state manipulator for Matlab
        template <typename ProblemClass>
        struct StateManipulator :
            public Optizelle::StateManipulator <ProblemClass>,
            public mxArrayPtr 
        {
        private:
            // Messaging object for reporting errors
            Messaging msg;

            // Keep a copy of a Matlab state lying around so that we can
            // use it to pass information back and forth to the Matlab
            // StateManipulator.
            mutable State <ProblemClass> mxstate;

            // Similarly, we keep only the Matlab versin of the bundle of
            // functions lying around
            mutable mxArrayPtr mxfns;

        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(StateManipulator);

            // We need the Matlab state manipulator, a copy of a Matlab state
            // to pass information, and a copy of the Matlab functions.
            StateManipulator(
                mxArray * const msg_,
                mxArray * const mxstate_,
                mxArray * const mxfns_,
                mxArray * const smanip_,
                mxArrayPtrMode::t const mode = mxArrayPtrMode::Capture
            ) :
                mxArrayPtr(smanip_,mode),
                msg(msg_,mxArrayPtrMode::Attach),
                mxstate(mxstate_,mxArrayPtrMode::Attach),
                mxfns(mxfns_,mxArrayPtrMode::Attach)
            {}

            // Application
            void eval(
                const typename ProblemClass::Functions::t & fns,
                typename ProblemClass::State::t & state,
                OptimizationLocation::t const & loc_
            ) const {
                // Convert the C++ state to a Matlab state
                mxstate.toMatlab(state);

                // Convert the lcoation to Matlab
                mxArrayPtr loc(OptimizationLocation::toMatlab(loc_));

                // Call the Matlab state manipulator give it mxstate and mxfns. 
                mxArray * eval(mxGetField(ptr,0,"eval"));
                std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject3(
                    eval,
                    mxfns.get(),
                    mxstate.get(),
                    loc.get()));
            
                // Check errors
                if(ret_err.second!=0)
                    msg.error("Evaluation of the StateManipulator object "
                        "failed.");

                // Convert the returned state to the C++ state 
                mxstate.reset(ret_err.first.release());
                mxstate.fromMatlab(state);
            }
        };

        // Vector space that works with Matlab objects
        template <typename Real=double> 
        struct MatlabVS { 
            // Prevent constructors 
            NO_CONSTRUCTORS(MatlabVS);

            // Setup the vector 
            typedef Optizelle::Matlab::Vector Vector; 

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
        typedef Optizelle::ScalarValuedFunction <double,MatlabVS>
            MxScalarValuedFunction;
        typedef Optizelle::VectorValuedFunction <double,MatlabVS,MatlabVS>
            MxVectorValuedFunction;
        typedef Optizelle::Operator <double,MatlabVS,MatlabVS> MxOperator;

        typedef Optizelle::Unconstrained <double,MatlabVS>
            MxUnconstrained;
        typedef Optizelle::EqualityConstrained <double,MatlabVS,MatlabVS>
            MxEqualityConstrained;
        typedef Optizelle::InequalityConstrained <double,MatlabVS,MatlabVS>
            MxInequalityConstrained;
        typedef Optizelle::Constrained <double,MatlabVS,MatlabVS,MatlabVS>
            MxConstrained;
        
        typedef Optizelle::json::Unconstrained <double,MatlabVS>
            MxJsonUnconstrained;
        typedef Optizelle::json::EqualityConstrained <double,MatlabVS,MatlabVS>
            MxJsonEqualityConstrained;
        typedef Optizelle::json::InequalityConstrained<double,MatlabVS,MatlabVS>
            MxJsonInequalityConstrained;
        typedef Optizelle::json::Constrained <double,MatlabVS,MatlabVS,MatlabVS>
            MxJsonConstrained;

        typedef typename Optizelle::RestartPackage <double>::t Reals;
        typedef typename Optizelle::RestartPackage <Natural>::t Naturals;
        typedef typename Optizelle::RestartPackage <std::string>::t Params;
        typedef typename Optizelle::RestartPackage <Vector>::t Vectors;

        // A simple scalar valued function interface, f : X -> R
        struct ScalarValuedFunction :
            public Optizelle::ScalarValuedFunction <double,MatlabVS>,
            public mxArrayPtr
        {
        private:
            // Create some type shortcuts
            typedef MatlabVS <>::Vector Vector; 

            // Messaging object
            Messaging msg; 

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(ScalarValuedFunction);

            // Create a function 
            explicit ScalarValuedFunction(
                mxArray * const msg_,
                mxArray * const f,
                mxArrayPtrMode::t mode=mxArrayPtrMode::Capture);

            // <- f(x) 
            double eval(Vector const & x) const; 

            // g = grad f(x) 
            void grad(Vector const & x,Vector & g) const; 

            // H_dx = hess f(x) dx 
            void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)const;
        };

        // A simple vector valued function interface, f : X -> Y
        struct VectorValuedFunction :
            public Optizelle::VectorValuedFunction<double,MatlabVS,MatlabVS>,
            public mxArrayPtr
        {
        private:
            // Create some type shortcuts
            typedef MatlabVS <>::Vector X_Vector; 
            typedef MatlabVS <>::Vector Y_Vector; 

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
                mxArray * const msg_,
                mxArray * const f,
                mxArrayPtrMode::t mode=mxArrayPtrMode::Capture);

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
            public Optizelle::Operator <double,MatlabVS,MatlabVS>,
            public mxArrayPtr
        {
        private:
            // Create some type shortcuts
            typedef MatlabVS <>::Vector X_Vector; 
            typedef MatlabVS <>::Vector Y_Vector; 

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
            mutable State <ProblemClass> mxstate;
            typename ProblemClass::State::t const & state;

            // Name of this function
            std::string const name;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Operator);

            // Create an operator 
            explicit Operator(
                std::string const & name_,
                mxArray * const msg_,
                mxArray * const op,
                mxArray * const mxstate_,
                typename ProblemClass::State::t const & state_,
                mxArrayPtrMode::t mode=mxArrayPtrMode::Capture
            ) :
                mxArrayPtr(op,mode),
                msg(msg_,mxArrayPtrMode::Attach),
                mxstate(mxstate_,mxArrayPtrMode::Attach),
                state(state_),
                name(name_)
            {}

            // y = A(x)
            void eval(X_Vector const & x,Y_Vector & y) const {
                // Convert the state to a Matlab state
                mxstate.toMatlab(state);

                // Apply the operator to the state, x, and y
                mxArray * eval(mxGetField(ptr,0,"eval"));
                std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                    eval,
                    mxstate.get(),
                    const_cast <X_Vector &> (x).get()));
                
                // Check errors
                if(ret_err.second!=0) {
                    std::stringstream ss;
                    ss << "Evaluation of the eval function in the operator "
                        << name << " failed.";
                    msg.error(ss.str());
                }
            
                // Copy the returned value into y 
                y.reset(ret_err.first.release());
            }
        };
        
        // Converts elements from C++ to Matlab 
        namespace toMatlab {
            // Sets a real in a Matlab state 
            void Real(
                std::string const & name,
                double const & value,
                mxArray * const obj 
            );

            // Sets a natural in a Matlab state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                mxArray * const obj 
            );
           
            // Sets a parameter in a Matlab state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<mxArray *(enum_t const &)> const & toMatlab,
                enum_t const & value,
                mxArray * const obj
            ) {
                mxArrayPtr item(toMatlab(value));
                mxSetField(obj,0,name.c_str(),item.release());
            }
            
            // Sets a vector in a Matlab state 
            void Vector(
                std::string const & name,
                Matlab::Vector const & value,
                mxArray * const obj 
            );
            
            // Sets a list of vectors in a Matlab state 
            void VectorList(
                std::string const & name,
                std::list <Matlab::Vector> const & values,
                mxArray * const obj 
            );
        
            // Sets restart vectors in Matlab 
            void Vectors(
                Matlab::Vectors const & values,
                mxArray * const mxvalues 
            );

            // Sets restart reals in Matlab 
            void Reals(
                Matlab::Reals const & values,
                mxArray * const mxvalues 
            );
            
            // Sets restart naturals in Matlab 
            void Naturals(
                Matlab::Naturals const & values,
                mxArray * const mxvalues 
            );
            
            // Sets restart parameters in Matlab 
            void Params(
                Matlab::Params const & values,
                mxArray * const mxvalues 
            );
        }

        // Converts elements from Matlab to C++ 
        namespace fromMatlab {
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                mxArray * const obj,
                double & value
            );
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                mxArray * const obj,
                Optizelle::Natural & value
            );
           
            // Sets a param C++ state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<enum_t(mxArray * const)> const & fromMatlab,
                mxArray * const obj,
                enum_t & value
            ) {
                mxArray * item(mxGetField(
                    const_cast <mxArray *> (obj),0,name.c_str()));
                value = fromMatlab(item);
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                mxArray * const obj,
                Matlab::Vector & value
            );
            
            // Sets a list of vectors in a C++ state 
            void VectorList(
                std::string const & name,
                mxArray * const obj,
                Matlab::Vector const & vec,
                std::list <Matlab::Vector> & values
            );
            
            // Sets a scalar-valued function in a C++ function bundle 
            void ScalarValuedFunction(
                std::string const & name,
                mxArray * const msg,
                mxArray * const obj,
                std::unique_ptr <MxScalarValuedFunction> & value
            );
            
            // Sets a vector-valued function in a C++ function bundle 
            void VectorValuedFunction(
                std::string const & name,
                mxArray * const msg,
                mxArray * const obj,
                std::unique_ptr <MxVectorValuedFunction> & value
            );
            
            // Sets an operator in a C++ function bundle 
            template <typename ProblemClass>
            void Operator(
                std::string const & name,
                mxArray * const msg,
                mxArray * const obj,
                mxArray * const mxstate,
                typename ProblemClass::State::t const & state,
                std::unique_ptr <MxOperator> & value
            ) {
                value.reset(new Matlab::Operator <ProblemClass> (
                    name,
                    msg,
                    mxGetField(obj,0,name.c_str()),
                    mxstate,
                    state,
                    mxArrayPtrMode::Attach));
            }
        
            // Sets restart vectors in C++ 
            void Vectors(
                Matlab::Vector const & vec,
                mxArray * const mxvalues,
                Matlab::Vectors & values
            );
            
            // Sets restart reals in C++ 
            void Reals(
                mxArray * const mxvalues,
                Matlab::Reals & values
            );
            
            // Sets restart naturals in C++ 
            void Naturals(
                mxArray * const mxvalues,
                Matlab::Naturals & values
            );
            
            // Sets restart parameters in C++ 
            void Params(
                mxArray * const mxvalues,
                Matlab::Params & values
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_();
                std::vector <char const *> fieldNames();
                
                // Create the structure for a Matlab state
                mxArray * mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxUnconstrained::State::t const & state,
                    mxArray * const mxstate
                );
                void toMatlab(
                    typename MxUnconstrained::State::t const & state,
                    mxArray * const mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArray * const mxstate,
                    typename MxUnconstrained::State::t & state
                );
                void fromMatlab(
                    mxArray * const mxstate,
                    typename MxUnconstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
                
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                template <typename ProblemClass>
                void fromMatlab_(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxUnconstrained::Functions::t & fns 
                ) {
                    fromMatlab::ScalarValuedFunction("f",msg,mxfns,fns.f);
                    fromMatlab::Operator <ProblemClass> (
                        "PH",msg,mxfns,mxstate,state,fns.PH);
                }
                void fromMatlab(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename MxUnconstrained::State::t const & state,
                    typename MxUnconstrained::Functions::t & fns 
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }
        
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_();
                std::vector <char const *> fieldNames();
                
                // Create the structure for a Matlab state
                mxArray * mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArray * const mxstate
                );
                void toMatlab(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArray * const mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArray * const mxstate,
                    typename MxEqualityConstrained::State::t & state
                );
                void fromMatlab(
                    mxArray * const mxstate,
                    typename MxEqualityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
                
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                template <typename ProblemClass>
                void fromMatlab_(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxEqualityConstrained::Functions::t & fns 
                ) {
                    fromMatlab::VectorValuedFunction("g",msg,mxfns,fns.g);
                    fromMatlab::Operator <ProblemClass> ("PSchur_left",
                        msg,mxfns,mxstate,state,fns.PSchur_left);
                    fromMatlab::Operator <ProblemClass> ("PSchur_right",
                        msg,mxfns,mxstate,state,fns.PSchur_right);
                }
                void fromMatlab(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename MxEqualityConstrained::State::t const & state,
                    typename MxEqualityConstrained::Functions::t & fns 
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }
        
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames_();
                std::vector <char const *> fieldNames();
                
                // Create the structure for a Matlab state
                mxArray * mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArray * const mxstate
                );
                void toMatlab(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArray * const mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArray * const mxstate,
                    typename MxInequalityConstrained::State::t & state
                );
                void fromMatlab(
                    mxArray * const mxstate,
                    typename MxInequalityConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
                
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                template <typename ProblemClass>
                void fromMatlab_(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxInequalityConstrained::Functions::t & fns 
                ) {
                    fromMatlab::VectorValuedFunction("h",msg,mxfns,fns.h);
                }
                void fromMatlab(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename MxInequalityConstrained::State::t const & state,
                    typename MxInequalityConstrained::Functions::t & fns 
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }
        
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
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
                // Returns the fields names of the state
                std::vector <char const *> fieldNames();
                
                // Create the structure for a Matlab state
                mxArray * mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab(
                    typename MxConstrained::State::t const & state,
                    mxArray * const mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab(
                    mxArray * const mxstate,
                    typename MxConstrained::State::t & state
                );
                
                // Creates a state and inserts the elements into mxstate 
                void create(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
                
                // Read json parameters from file
                void readJson(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }

            // All the functions required by an optimization algorithm.  
            namespace Functions{
                // Convert a Matlab bundle to C++ 
                void fromMatlab(
                    mxArray * const msg,
                    mxArray * const mxfns,
                    mxArray * const mxstate,
                    typename MxConstrained::State::t const & state,
                    typename MxConstrained::Functions::t & fns 
                );
            }

            // This contains the different algorithms used for optimization
            namespace Algorithms {
                // Solves an optimization problem
                void getMin(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }
        
            // Utilities for restarting the optimization
            namespace Restart {
                // Release the data into structures controlled by the user 
                void release(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Capture data from structures controlled by the user.  
                void capture(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Writes a json restart file
                void write_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );

                // Reads a json restart file
                void read_restart(
                    int nOutput,mxArray* pOutput[],
                    int nInput,mxArray* pInput[]
                );
            }
        }
    }
}

#endif
