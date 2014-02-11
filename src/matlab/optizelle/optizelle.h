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
            
                // Call the Matlab state manipulator 
                // give it mxstate and mxfns.  Note, mxfns is given raw.
                mxArrayPtr eval(mxGetField(ptr,0,"eval"));
                std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject3(
                    eval.get(),
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
                mxArrayPtr eval(mxGetField(ptr,0,"eval"));
                std::pair <mxArrayPtr,int> ret_err(mxArray_CallObject2(
                    eval.get(),
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
                mxArrayPtr item(mxGetField(
                    const_cast <mxArray *> (obj),0,name.c_str()));
                value = fromMatlab(item.get());
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
                    state));
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
    }
}

#if 0
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
void handleException(int err,std::string msg) {
    if(err) {
        // Print out our error message and then return control
        // to Matlab.
        std::stringstream ss;
        ss << msg << std::endl; 
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
    mxArrayPtr rand_;
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
        axpy_(NULL), innr_(NULL), rand_(NULL), prod_(NULL), id_(NULL),
        linv_(NULL), barr_(NULL), srch_(NULL), symm_(NULL) {}

    // On a simple vector, initialize both the internal storage as well as
    // the basic linear algebra.
    MatVector(
        mxArray* vec,
        mxArray* copy__,
        mxArray* scal__,
        mxArray* zero__,
        mxArray* axpy__,
        mxArray* innr__,
        mxArray* rand__ 
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
        rand_.reset(mxDuplicateArray(rand__));
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
        mxArray* rand__,
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
        rand_.reset(mxDuplicateArray(rand__));
        prod_.reset(mxDuplicateArray(prod__));
        id_.reset(mxDuplicateArray(id__));
        linv_.reset(mxDuplicateArray(linv__));
        barr_.reset(mxDuplicateArray(barr__));
        srch_.reset(mxDuplicateArray(srch__));
        symm_.reset(mxDuplicateArray(symm__));
    }

    // Create a move constructor so we can interact with stl objects
    MatVector(MatVector && vec) noexcept {
        // Grab the memory for the vector 
        reset(vec.release());

        // Now, grab all of the algebra functions
        copy_.reset(vec.copy_.release());
        scal_.reset(vec.scal_.release());
        zero_.reset(vec.zero_.release());
        axpy_.reset(vec.axpy_.release());
        innr_.reset(vec.innr_.release());
        rand_.reset(vec.rand_.release());
        prod_.reset(vec.prod_.release());
        id_.reset(vec.id_.release());
        linv_.reset(vec.linv_.release());
        barr_.reset(vec.barr_.release());
        srch_.reset(vec.srch_.release());
        symm_.reset(vec.symm_.release());
    }
    
    // Move assignment operator
    MatVector const & operator = (MatVector && vec) noexcept {
        // Grab the memory for the vector 
        reset(vec.release());

        // Now, grab all of the algebra functions
        copy_.reset(vec.copy_.release());
        scal_.reset(vec.scal_.release());
        zero_.reset(vec.zero_.release());
        axpy_.reset(vec.axpy_.release());
        innr_.reset(vec.innr_.release());
        rand_.reset(vec.rand_.release());
        prod_.reset(vec.prod_.release());
        id_.reset(vec.id_.release());
        linv_.reset(vec.linv_.release());
        barr_.reset(vec.barr_.release());
        srch_.reset(vec.srch_.release());
        symm_.reset(vec.symm_.release());

        return *this;
    }

    // Memory allocation and size setting 
    MatVector init() {
        // Create and return a vector based on the internal
        if( prod_.get() && id_.get() && linv_.get() && barr_.get() &&
            srch_.get() && symm_.get()
        )
            return std::move(MatVector(
                get(),
                copy_.get(),
                scal_.get(),
                zero_.get(),
                axpy_.get(),
                innr_.get(),
                rand_.get(),
                prod_.get(),
                id_.get(),
                linv_.get(),
                barr_.get(),
                srch_.get(),
                symm_.get()));
        else
            return std::move(MatVector(
                get(),
                copy_.get(),
                scal_.get(),
                zero_.get(),
                axpy_.get(),
                innr_.get(),
                rand_.get()));
    }

    // y <- x (Shallow.  No memory allocation.) 
    void copy(const mxArray* x){
        reset(mxDuplicateArray(x));
    }

    // y <- x (Shallow.  No memory allocation.) 
    void copy(const MatVector& x) {
        mxArray* input[2]={copy_.get(),x.get()};
        mxArray* output[1];
        int err=mexCallMATLAB(1,output,2,input,"feval");
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
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,"Error in the vector space scal function");
        reset(output[0]);
    }

    // x <- 0 
    void zero() {
        // Create the inputs and outputs 
        mxArray* input[2]={zero_.get(),get()};
        mxArray* output[1];

        // Find the zero vector and store the result.
        int err=mexCallMATLAB(1,output,2,input,"feval");
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
        int err=mexCallMATLAB(1,output,4,input,"feval");
        handleException(err,"Error in the vector space axpy function");
        reset(output[0]);
    }

    // innr <- <x,y> 
    double innr(const MatVector& x) const {
        // Create the inputs and outputs 
        mxArray* input[3]={innr_.get(),x.get(),get()};
        mxArray* output[1];

        // Compute the inner product 
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,"Error in the vector space innr function");
        mxArrayPtr output_(output[0]);

        // Get the result of the computation 
        double alpha=mxGetScalar(output_.get());

        // Return the result 
        return alpha;
    }

    // x <- random 
    void rand() {
        // Create the inputs and outputs 
        mxArray* input[2]={rand_.get(),get()};
        mxArray* output[1];

        // Find the rand vector and store the result.
        int err=mexCallMATLAB(1,output,2,input,"feval");
        handleException(err,"Error in the vector space rand function");
        reset(output[0]);
    }

    // Jordan product, z <- x o y 
    void prod(const MatVector& x, const MatVector& y) {
        // Create the inputs and outputs 
        mxArray* input[3]={prod_.get(),x.get(),y.get()};
        mxArray* output[1];

        // Compute the product and store the result 
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,"Error in the vector space prod function");
        reset(output[0]);
    }

    // Identity element, x <- e such that x o e = x 
    void id() {
        // Create the inputs and outputs 
        mxArray* input[2]={id_.get(),get()};
        mxArray* output[1];

        // Find the identity element and store the result
        int err=mexCallMATLAB(1,output,2,input,"feval");
        handleException(err,"Error in the vector space id function");
        reset(output[0]);
    }

    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y 
    void linv(const MatVector& x, const MatVector& y) {
        // Create the inputs and outputs 
        mxArray* input[3]={linv_.get(),x.get(),y.get()};
        mxArray* output[1];

        // Compute the product inverse and store the result 
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,"Error in the vector space linv function");
        reset(output[0]);
    }

    // Barrier function, barr <- barr(x) where x o grad barr(x) = e 
    double barr() const {
        // Create the inputs and outputs 
        mxArray* input[2]={barr_.get(),get()};
        mxArray* output[1];

        // Compute the barrier function 
        int err=mexCallMATLAB(1,output,2,input,"feval");
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
        int err=mexCallMATLAB(1,output,3,input,"feval");
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
        int err=mexCallMATLAB(1,output,2,input,"feval");
        handleException(err,"Error in the vector space symm function");
        reset(output[0]);
    }
};

        
template <typename Real=double> 
struct MatlabVS { 
    // Setup the vector 
    typedef MatVector Vector; 
        
    // Memory allocation and size setting 
    static Vector init(const Vector & x) {
        return std::move(const_cast <Vector &> (x).init());
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
        
    // x <- random 
    static void rand(Vector& x) { 
        x.rand();
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
    : public Optizelle::ScalarValuedFunction <double,MatlabVS>
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
    double eval(const MatVector& x) const {
        // Create the inputs and outputs 
        mxArray* input[2]={eval_.get(),x.get()}; 
        mxArray* output[1];

        // Evaluate the function 
        mexCallMATLAB(1,output,2,input,"feval"); 
        int err=mexCallMATLAB(1,output,2,input,"feval");
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
        int err=mexCallMATLAB(1,output,2,input,"feval");
        handleException(err,
            "Error in the scalar valued function (objective) gradient");
        g.reset(output[0]);
    }

    // H_dx = hess f(x) dx 
    void hessvec(const MatVector& x,const MatVector& dx,MatVector& H_dx) const {
        mxArray* input[3]={hessvec_.get(),x.get(),dx.get()}; 
        mxArray* output[1];
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,
            "Error in the scalar valued function (objective) Hessian-vector "
            "product");
        H_dx.reset(output[0]);
    }
};

// A vector valued function using Matlab provided functions 
struct MatlabVectorValuedFunction
    : public Optizelle::VectorValuedFunction <double,MatlabVS,MatlabVS> {
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
    void eval(const X_Vector& x,Y_Vector& y) const { 
        mxArray* input[2]={eval_.get(),x.get()}; 
        mxArray* output[1];
        mexCallMATLAB(1,output,2,input,"feval");
        int err=mexCallMATLAB(1,output,2,input,"feval");
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
        int err=mexCallMATLAB(1,output,3,input,"feval");
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
        int err=mexCallMATLAB(1,output,3,input,"feval");
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
        int err=mexCallMATLAB(1,output,4,input,"feval");
        handleException(err,
            "Error in the vector valued function (constraint) second "
            "derivative adjoint");
        z.reset(output[0]);
     }
};

// An operator using Matlab provided functions 
struct MatlabOperator : public Optizelle::Operator <double,MatlabVS,MatlabVS> {
private:
    // Create some type shortcuts
    typedef MatlabVS <>::Vector X_Vector;
    typedef MatlabVS <>::Vector Y_Vector;

    // The Matlab functions
    mxArrayPtr eval_;

    // The base point around where we evaluate the operator
    X_Vector& x; 

    // Name of the operator
    std::string name;
    
public:
    // When we construct the function, we read in all the names from the
    // Matlab structure.  We need to check the structure prior to
    // constructing this object.
    MatlabOperator(mxArray* f,X_Vector& x_,std::string name_)
        :  eval_(nullptr), x(x_), name(name_)
    {
        eval_.reset(mxDuplicateArray(mxGetField(f,0,"eval"))); 
    }

    // y=f(x)
    virtual void eval(const X_Vector& dx,Y_Vector& y) const { 
        mxArray* input[3]={eval_.get(),x.get(),dx.get()}; 
        mxArray* output[1];
        int err=mexCallMATLAB(1,output,3,input,"feval");
        handleException(err,"Error in the operator " + name);
        y.reset(output[0]);
    }
};

// A messaging utility that hooks directly into Matlab
struct MatlabMessaging : public Optizelle::Messaging {
    // Prints a message
    void print(const std::string& msg) const {
        std::string msg_ = msg + "\n";
        mexPrintf(msg_.c_str());
        mexEvalString("drawnow");
    }

    // Prints an error
    void error(const std::string& msg) const {
        mexErrMsgTxt(msg.c_str());
    }
};

// Checks that a structure array contains certain functions
void check_fns(
    const mxArray* X,
    std::string ops[],
    Optizelle::Natural nops,
    std::string type,
    std::string name
) {
    // Create a type shortcut
    typedef Optizelle::Natural Natural;

    // The argument should be a structure array
    if(!mxIsStruct(X))
        MatlabMessaging().error("The " + type + " " + name +
            " must be a structure array.");

    // Check that the structure contains all the functions
    for(Natural i=0;i<nops;i++)
        if(mxGetField(X,0,ops[i].c_str())==NULL)
            MatlabMessaging().error("Missing the " + ops[i] + 
                " function in the " + type + " " + name +".");
    
    // Check that all of the listed functions are really functions 
    for(Natural i=0;i<nops;i++)
        if(!mxIsClass(mxGetField(X,0,ops[i].c_str()),"function_handle"))
            MatlabMessaging().error("The field " + ops[i] +
                " in the " + type + " " + name + " must be a function.");
}

// Check that we have all the operations necessary for a vector space
void check_vector_space(const mxArray* X,std::string name) {
    std::string ops[6] = { "copy", "scal", "zero", "axpy", "innr", "rand" };
    check_fns(X,ops,6,"vector space",name);
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
Optizelle::ProblemClass::t get_problem_class_from_vs(const mxArray* VS) {
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
    Optizelle::ProblemClass::t problem_class = Optizelle::ProblemClass::Unconstrained;

    // Next, check if we have any equality constraints
    if(mxGetField(VS,0,"Y") != NULL) {

        // Check for the vector space operations
        check_vector_space(mxGetField(VS,0,"Y"),"Y");

        // Adjust the algorithm class
        problem_class
            = Optizelle::ProblemClass::EqualityConstrained;
    }

    // Next, check if we have inequality constraints
    if(mxGetField(VS,0,"Z") != NULL) {

        // Check for the Euclidean-Jordan algebra operations
        check_euclidean_jordan(mxGetField(VS,0,"Z"),"Z");

        // Adjust the algorithm class
        if(problem_class==Optizelle::ProblemClass::EqualityConstrained)
            problem_class = Optizelle::ProblemClass::Constrained;
        else
            problem_class = Optizelle::ProblemClass::InequalityConstrained;
    }

    return problem_class;
}

// Determine the problem class from a bundle of vector spaces 
Optizelle::ProblemClass::t get_problem_class_from_fns(const mxArray* fns) {
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
    Optizelle::ProblemClass::t problem_class
        = Optizelle::ProblemClass::Unconstrained;

    // Next, check if we have any equality constraints
    if(mxGetField(fns,0,"g") != NULL) {

        // Check for the vector space operations
        check_vector_valued_fn(mxGetField(fns,0,"g"),"g");

        // Adjust the algorithm class
        problem_class
            = Optizelle::ProblemClass::EqualityConstrained;
    }

    // Next, check if we have inequality constraints
    if(mxGetField(fns,0,"h") != NULL) {

        // Check for the Euclidean-Jordan algebra operations
        check_vector_valued_fn(mxGetField(fns,0,"h"),"h");

        // Adjust the algorithm class
        if(problem_class==Optizelle::ProblemClass::EqualityConstrained)
            problem_class = Optizelle::ProblemClass::Constrained;
        else
            problem_class = Optizelle::ProblemClass::InequalityConstrained;
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
    Optizelle::ProblemClass::t problem_class_vs=get_problem_class_from_vs(VS);
    Optizelle::ProblemClass::t problem_class_fn=get_problem_class_from_fns(fns);

    // Make sure the problem class inferred by the vector space matches
    // that from the function class
    if(problem_class_vs != problem_class_fn) {
        std::stringstream ss;
        ss << "The problem class inferred from the vector space ("
            << Optizelle::ProblemClass::to_string(problem_class_vs) 
            << ") does not match\nthe problem class inferred from the "
            << "optimization functions ("
            << Optizelle::ProblemClass::to_string(problem_class_fn)
            << ").";

        MatlabMessaging().error(ss.str());
    }
}

// Check that the bundle of points aligns with the current vector space.
// We use this bundle of points for optimization 
void check_optimization_pts(const mxArray* VS,const mxArray* pts) {
    // Get the problem class from the vector space
    Optizelle::ProblemClass::t problem_class=get_problem_class_from_vs(VS);
    
    // First, check that the bundle of points is a structure array        
    if(!mxIsStruct(pts))
        MatlabMessaging().error("The bundle of points must be a "
            "structure array.");

    // Next, insure that we have x, y, and z 
    if(mxGetField(pts,0,"x")==NULL)
        MatlabMessaging().error("The bundle of points must contain x.");
    if((problem_class == Optizelle::ProblemClass::EqualityConstrained 
        || problem_class == Optizelle::ProblemClass::Constrained) 
        && mxGetField(pts,0,"y")==NULL
    )
        MatlabMessaging().error("The bundle of points must contain y.");
    if((problem_class == Optizelle::ProblemClass::InequalityConstrained 
        || problem_class == Optizelle::ProblemClass::Constrained) 
        && mxGetField(pts,0,"z")==NULL
    )
        MatlabMessaging().error("The bundle of points must contain z.");
}

// Check that the bundle of points aligns with the current vector space.
// We use this bundle of points for the finite difference tests.
void check_diagnostic_pts(const mxArray* VS,const mxArray* pts) {
    // Get the problem class from the vector space
    Optizelle::ProblemClass::t problem_class=get_problem_class_from_vs(VS);

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
    if(problem_class == Optizelle::ProblemClass::EqualityConstrained ||
       problem_class == Optizelle::ProblemClass::Constrained
    ) {
        if(mxGetField(pts,0,"dy")==NULL)
            MatlabMessaging().error("For equality constrained problems, "
                "the bundle of points must contain dy.");
    }

    // If we're inequality constrained, check for dz
    if(problem_class == Optizelle::ProblemClass::InequalityConstrained ||
       problem_class == Optizelle::ProblemClass::Constrained
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
        MatlabMessaging().error("Optizelle usage:\n"
            "Optimization: x=Optizelle(VS,fns,pts,params)\n"
            "Diagnostics:  Optizelle(VS,fns,pts)");

    // Check the vector spaces
    checkVS(pInput[0]);

    // Check the bundle of functions
    checkFns(pInput[0],pInput[1]);

    // Determine if we're optimizing or doing diagnostics
    bool diagnostics = nInput==3;
     
    // Get the problem class
    Optizelle::ProblemClass::t problem_class=get_problem_class_from_vs(pInput[0]);

    // Create a bundle of functions
    Optizelle::Constrained <double,MatlabVS,MatlabVS,MatlabVS>::Functions::t fns;
    fns.f.reset(new MatlabScalarValuedFunction (mxGetField(pInput[1],0,"f")));
    switch(problem_class) {
    case Optizelle::ProblemClass::Unconstrained:
        break;
    case Optizelle::ProblemClass::EqualityConstrained:
        fns.g.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"g")));
        break;
    case Optizelle::ProblemClass::InequalityConstrained:
        fns.h.reset(new MatlabVectorValuedFunction
            (mxGetField(pInput[1],0,"h")));
        break;
    case Optizelle::ProblemClass::Constrained:
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
    std::unique_ptr <MatVector> x;
    std::unique_ptr <MatVector> dx;
    std::unique_ptr <MatVector> dxx;
    std::unique_ptr <MatVector> y;
    std::unique_ptr <MatVector> dy;
    std::unique_ptr <MatVector> z;
    std::unique_ptr <MatVector> dz;

    // Check that we have a valid set of optimization points (x,y,z)
    check_optimization_pts(pInput[0],pInput[2]);

    // Always initialize x
    x.reset(new MatVector(
        mxGetField(pInput[2],0,"x"),
        mxGetField(X,0,"copy"),
        mxGetField(X,0,"scal"),
        mxGetField(X,0,"zero"),
        mxGetField(X,0,"axpy"),
        mxGetField(X,0,"innr"),
        mxGetField(X,0,"rand")));
    
    // If we need to do a finite difference, initialize dx and dxx
    if(diagnostics) {
        // Check that we have all the valid directions for the finite
        // difference tests.
        check_diagnostic_pts(pInput[0],pInput[2]);

        // Get the points for checking the objective
        dx.reset(new MatVector(x->init()));
            dx->copy(mxGetField(pInput[2],0,"dx"));
        dxx.reset(new MatVector(x->init()));
            dxx->copy(mxGetField(pInput[2],0,"dx"));
    }
    
    // Possibly initialize y and dy
    if(problem_class == Optizelle::ProblemClass::EqualityConstrained ||
       problem_class == Optizelle::ProblemClass::Constrained
    ) {
        mxArray* Y(mxGetField(pInput[0],0,"Y"));
        y.reset(new MatVector(
            mxGetField(pInput[2],0,"y"),
            mxGetField(Y,0,"copy"),
            mxGetField(Y,0,"scal"),
            mxGetField(Y,0,"zero"),
            mxGetField(Y,0,"axpy"),
            mxGetField(Y,0,"innr"),
            mxGetField(Y,0,"rand")));
        if(diagnostics) {
            dy.reset(new MatVector(y->init()));
                dy->copy(mxGetField(pInput[2],0,"dy"));
        }
    }

    // Possibly initialize z and dz
    if(problem_class == Optizelle::ProblemClass::InequalityConstrained ||
       problem_class == Optizelle::ProblemClass::Constrained
    ) {
        mxArray* Z(mxGetField(pInput[0],0,"Z"));
        z.reset(new MatVector(
            mxGetField(pInput[2],0,"z"),
            mxGetField(Z,0,"copy"),
            mxGetField(Z,0,"scal"),
            mxGetField(Z,0,"zero"),
            mxGetField(Z,0,"axpy"),
            mxGetField(Z,0,"innr"),
            mxGetField(Z,0,"rand"),
            mxGetField(Z,0,"prod"),
            mxGetField(Z,0,"id"),
            mxGetField(Z,0,"linv"),
            mxGetField(Z,0,"barr"),
            mxGetField(Z,0,"srch"),
            mxGetField(Z,0,"symm")));
        if(diagnostics) {
            dz.reset(new MatVector(z->init()));
                dz->copy(mxGetField(pInput[2],0,"dz"));
        }
    }

    // If we're doing diagnostics, diagnose
    if(diagnostics) {
        // Do a finite difference check and symmetric check on the objective
        MatlabMessaging().print("Diagnostics on the objective.");
        Optizelle::Diagnostics::gradientCheck<double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx,"f");
        Optizelle::Diagnostics::hessianCheck <double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx,"f");
        Optizelle::Diagnostics::hessianSymmetryCheck <double,MatlabVS>
            (MatlabMessaging(),*(fns.f),*x,*dx,*dxx,"f");

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == Optizelle::ProblemClass::EqualityConstrained ||
           problem_class == Optizelle::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the equality constraint.");
        
            // Do some finite difference tests on the constraint
            Optizelle::Diagnostics::derivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy,"g");
            Optizelle::Diagnostics
                ::derivativeAdjointCheck<double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy,"g");
            Optizelle::Diagnostics
                ::secondDerivativeCheck<double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.g),*x,*dx,*dy,"g");
        }

        // Run diagnostics on the equality constraints if necessary 
        if(problem_class == Optizelle::ProblemClass::InequalityConstrained ||
           problem_class == Optizelle::ProblemClass::Constrained
        ) {
            MatlabMessaging()
                .print("\nDiagnostics on the inequality constraint.");

            // Do some finite difference tests on the constraint
            Optizelle::Diagnostics::derivativeCheck <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz,"h");
            Optizelle::Diagnostics::derivativeAdjointCheck
                <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz,"h");
            Optizelle::Diagnostics::secondDerivativeCheck
                <double,MatlabVS,MatlabVS>
                (MatlabMessaging(),*(fns.h),*x,*dx,*dz,"h");
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
        case Optizelle::ProblemClass::Unconstrained: {
            // Allocate memory for the state
            Optizelle::Unconstrained <double,MatlabVS>::State::t state(*x);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x,"PH"));

            // Read the parameters and optimize
            Optizelle::json::Unconstrained <double,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            Optizelle::Unconstrained<double,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);

            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                Optizelle::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case Optizelle::ProblemClass::InequalityConstrained: {
            // Allocate memory for the state
            Optizelle::InequalityConstrained <double,MatlabVS,MatlabVS>::State::t
                state(*x,*z);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x,"PH"));

            // Read the parameters and optimize
            Optizelle::json::InequalityConstrained <double,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            Optizelle::InequalityConstrained<double,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.release());
            mxSetField(pOutput[0],0,"z",state.z.release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                Optizelle::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case Optizelle::ProblemClass::EqualityConstrained: {
            // Allocate memory for the state
            Optizelle::EqualityConstrained <double,MatlabVS,MatlabVS>::State::t
                state(*x,*y);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x,"PH"));
            if(mxGetField(pInput[1],0,"PSchur_left")!=NULL)
                fns.PSchur_left.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_left"),state.x,
                    "PSchur_left"));
            if(mxGetField(pInput[1],0,"PSchur_right")!=NULL)
                fns.PSchur_right.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_right"),state.x,
                    "PSchur_right"));

            // Read the parameters and optimize
            Optizelle::json::EqualityConstrained <double,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            Optizelle::EqualityConstrained<double,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.release());
            mxSetField(pOutput[0],0,"y",state.y.release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                Optizelle::StoppingCondition::to_string(state.opt_stop).c_str()));
            break;
        } case Optizelle::ProblemClass::Constrained: {
            // Allocate memory for the state
            Optizelle::Constrained <double,MatlabVS,MatlabVS,MatlabVS>::State::t
                state(*x,*y,*z);

            // If we have a preconditioner, add it
            if(mxGetField(pInput[1],0,"PH")!=NULL)
                fns.PH.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PH"),state.x,"PH"));
            if(mxGetField(pInput[1],0,"PSchur_left")!=NULL)
                fns.PSchur_left.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_left"),state.x,
                    "PSchur_left"));
            if(mxGetField(pInput[1],0,"PSchur_right")!=NULL)
                fns.PSchur_right.reset(new MatlabOperator
                    (mxGetField(pInput[1],0,"PSchur_right"),state.x,
                    "PSchur_right"));

            // Read the parameters and optimize
            Optizelle::json::Constrained <double,MatlabVS,MatlabVS,MatlabVS>
                ::read(MatlabMessaging(),params,state);
            Optizelle::Constrained<double,MatlabVS,MatlabVS,MatlabVS>::Algorithms
                ::getMin(MatlabMessaging(),fns,state);
            
            // Save the answer
            mxSetField(pOutput[0],0,"x",state.x.release());
            mxSetField(pOutput[0],0,"y",state.y.release());
            mxSetField(pOutput[0],0,"z",state.z.release());
            mxSetField(pOutput[0],0,"opt_stop",mxCreateString(
                Optizelle::StoppingCondition::to_string(state.opt_stop).c_str()));
        } }
    }
}
#endif
