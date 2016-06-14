/*
Copyright 2013-2014 OptimoJoe.

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
    // Forward declare some pointer types
    namespace Matlab {
        struct mxArrayPtr;
        struct mxManaged;
        struct mxUnmanaged;
    }

    // Extend our enumerated types to convert to and from MATLAB/Octave
    namespace OptimizationStop { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & opt_stop);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }

    namespace TruncatedStop {
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & trunc_stop);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace AlgorithmClass { 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & algorithm_class);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace Operators{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & op);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace LineSearchDirection{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & dir);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace LineSearchKind{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & kind);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace OptimizationLocation{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & loc);

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace FunctionDiagnostics{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & diag); 

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace VectorSpaceDiagnostics{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & diag); 

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    
    namespace DiagnosticScheme{ 
        // Converts t to a Matlab enumerated type
        Matlab::mxManaged toMatlab(t const & dscheme); 

        // Converts a Matlab enumerated type to t 
        t fromMatlab(Matlab::mxArrayPtr const & member);
    }
    namespace Matlab {
        // Exception for when a function in MATLAB/Octave throws an error
        namespace Exception {
            struct t : public std::runtime_error {
                using std::runtime_error::runtime_error;
            };
        }

        // Wrapper for mxArray *
        struct mxArrayPtr {
        protected:
            // Internal storage of the pointer
            mxArray * ptr;
            
            // Grab the pointer.  Keep this protected to prevent direct
            // instantiation of this class.
            mxArrayPtr(
                mxArray const * const & ptr_
            );

        public:
            // Allow move semantics
            mxArrayPtr(mxArrayPtr && p);
            mxArrayPtr & operator = (mxArrayPtr && p);

            // Disallow copy semantics
            mxArrayPtr(mxArrayPtr const &) = delete;
            mxArrayPtr & operator = (mxArrayPtr const &) = delete;

            // Allow the deallocator to be overwritten 
            virtual ~mxArrayPtr();

            // Grab the internal pointer
            mxArray * get() const;

            // Release ownership and management of the pointer
            mxArray * release();
        };

        // Unmanaged pointers.  This assumes that MATLAB/Octave will free the
        // memory.
        struct mxUnmanaged : public mxArrayPtr {
            // Grab the pointer
            mxUnmanaged(
                mxArray const * const & ptr
            ); 

            // Allow copy and move semantics
            mxUnmanaged(mxUnmanaged const & p);
            mxUnmanaged & operator = (mxUnmanaged const & p);
            mxUnmanaged(mxUnmanaged &&) = default;
            mxUnmanaged & operator = (mxUnmanaged &&) = default;
        };

        // Managed pointers.  We free the memory, not MATLAB/Octave.
        struct mxManaged : public mxArrayPtr {
            // Grab the pointer
            mxManaged(
                mxArray const * const & ptr
            );

            // Free the memory
            virtual ~mxManaged();

            // Allow move semantics
            mxManaged(mxManaged &&) = default;
            mxManaged & operator = (mxManaged &&) = default;

            // Reset the pointer to a given value
            void reset(mxArray * ptr_);
        };

        // Calls a Matlab function with one argument and no returns 
        int mxArray_CallObject1_0(
            mxArray * const fn,
            mxArray * const arg1
        );

        // Calls a Matlab function with one argument 
        std::tuple <mxManaged,int> mxArray_CallObject1(
            mxArray * const fn,
            mxArray * const arg1
        );
        
        // Calls a Matlab function with two arguments
        std::tuple <mxManaged,int> mxArray_CallObject2(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2
        );
        
        // Calls a Matlab function with three arguments
        std::tuple <mxManaged,int> mxArray_CallObject3(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3
        );
        
        // Calls a Matlab function with four arguments
        std::tuple <mxManaged,int> mxArray_CallObject4(
            mxArray * const fn,
            mxArray * const arg1,
            mxArray * const arg2,
            mxArray * const arg3,
            mxArray * const arg4
        );

        // Creates a Matlab double from a C++ double
        mxManaged mxArray_FromDouble(double const x_);

        // Creates a Matlab int from a C++ size_t 
        mxManaged mxArray_FromSize_t(Natural const x_);

        // Converts an Optizelle enumerated type to a mxArray *
        mxManaged enumToMxArray(
            std::string const & type,
            std::string const & member 
        );
        
        // Converts an Optizelle enumerated type to a Natural
        Natural enumToNatural(
            std::string const & type,
            std::string const & member 
        );

        // Converts a MATLAB double to an Optizelle Natural
        Natural fromDouble(double value);

        // A messaging utility that hooks directly into MATLAB/Octave 
        namespace Messaging {
            template <typename A>
            Optizelle::Messaging::t matlab(A && print_) {
                auto print = std::shared_ptr <mxArrayPtr> (
                    new A(std::move(print_)));
                return [print](std::string const & msg_) {

                    // Call the print function
                    auto msg = mxManaged(mxCreateString(msg_.c_str()));
                    auto err = mxArray_CallObject1_0(
                        print->get(),
                        msg.get());

                    // Check errors
                    if(err)
                        throw Matlab::Exception::t( __LOC__
                            + ", evaluation of the Messaging function failed");
                };
            }
        }

        // This class merges the vector space with a vector into a singular 
        // object.  We require this structure since Optizelle requires the
        // vector space to be static.  Since the user is passing us a vector
        // space dynamically, we merge the vector space functions with the
        // vectors and then statically define the vector space to call these
        // functions.
        struct Vector {
        private:
            // Vector space
            std::shared_ptr <mxArrayPtr> vs;

        public:
            // Data
            std::unique_ptr <mxArrayPtr> data;

            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Vector)

            // Build a vector

            // Generally, we use this version for the initial vector space
            template <typename A,typename B>
            Vector(
                A && vs_,
                B && data_
            ) :
                vs(new A(std::move(vs_))),
                data(new B(std::move(data_)))
            {}

            // Then, we use this version to spawn new vectors
            template <typename A>
            Vector(
                std::shared_ptr <mxArrayPtr> const & vs_,
                A && data_
            ) :
                vs(vs_),
                data(new A(std::move(data_)))
            {}

            // Allow move constructors
            Vector(Vector &&) = default;
            Vector & operator = (Vector &&) = default;

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
            
            // Converts (copies) a value into Matlab.  
            mxManaged toMatlab() const;
            
            // Converts (copies) a value from Matlab.  This assumes that the
            // vector space functions have already been properly assigned.
            void fromMatlab(mxArrayPtr const & ptr);
        };

        // Matlab state
        template <typename ProblemClass>
        struct State {
            // Data
            std::unique_ptr <mxArrayPtr> data;

            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(State)

            // On construction, we just grab the pointer to the state object
            template <typename A>
            State(A && data_) :
                data(new A(std::move(data_)))
            {}

            // Allow move constructors
            State(State &&) = default;
            State & operator = (State &&) = default;

            // Convert a C++ state to a Matlab state 
            void toMatlab(typename ProblemClass::State::t const & state);

            // Convert a Matlab state to C++ 
            void fromMatlab(typename ProblemClass::State::t & state);
        };
        
        // Matlab bundle of functions 
        template <typename ProblemClass>
        struct Functions {
        private:
            // Keep two states lying around, so that we can communicate this to
            // our operator.  Here, state is linked to the optimization state
            // and this is used to copy into mxstate, which is given to the
            // operator.
            State <ProblemClass> & mxstate;
            typename ProblemClass::State::t const & state;
            
        public:
            // Data
            std::unique_ptr <mxArrayPtr> data;

            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Functions)

            // On construction, we just grab the pointer to the bundle object
            template <typename A>
            Functions(
                State <ProblemClass> & mxstate_,
                typename ProblemClass::State::t const & state_,
                A && data_
            ) :
                mxstate(mxstate_),
                state(state_),
                data(new A(std::move(data_)))
            {}

            // Convert a Matlab bundle to C++ 
            void fromMatlab(typename ProblemClass::Functions::t & fns);
        };
        
        // The state manipulator for Matlab
        template <typename ProblemClass>
        struct StateManipulator :
            public Optizelle::StateManipulator <ProblemClass>
        {
        private:
            // Keep a copy of a Matlab state lying around so that we can
            // use it to pass information back and forth to the Matlab
            // StateManipulator
            State <ProblemClass> & mxstate;

            // Similarly, keep a copy of the functions lying around to pass to
            // the StateManipulator 
            Functions <ProblemClass> & mxfns;
            
            // Underlying state manipulator 
            std::unique_ptr <mxArrayPtr> smanip;

        public:
            // Disallow constructors
            NO_DEFAULT_COPY_ASSIGNMENT(StateManipulator)

            // We need the Matlab state manipulator, a copy of a Matlab state
            // to pass information, and a copy of the Matlab functions.
            template <typename A>
            StateManipulator(
                State <ProblemClass> & mxstate_,
                Functions <ProblemClass> & mxfns_,
                A && smanip_
            ) :
                mxstate(mxstate_),
                mxfns(mxfns_),
                smanip(new A(std::move(smanip_)))
            {}

            // Application
            void eval(
                typename ProblemClass::Functions::t const & fns,
                typename ProblemClass::State::t & state,
                OptimizationLocation::t const & loc_
            ) const {
                // Convert the C++ state to a Matlab state
                mxstate.toMatlab(state);

                // Convert the lcoation to Matlab
                auto loc = OptimizationLocation::toMatlab(loc_);

                // Call the Matlab state manipulator give it mxstate and mxfns. 
                auto eval = mxGetField(smanip->get(),0,"eval");
                auto ret_err = mxArray_CallObject3(
                    eval,
                    mxfns.data->get(),
                    mxstate.data->get(),
                    loc.get());
            
                // Check errors
                if(std::get <1> (ret_err))
                    throw Matlab::Exception::t( __LOC__
                        + ", evaluation of the StateManipulator object failed");

                // Convert the returned state to the C++ state 
                mxstate.data.reset(
                    new mxManaged(std::move(std::get<0> (ret_err))));
                mxstate.fromMatlab(state);
            }
        };

        // Vector space that works with Matlab objects
        template <typename Real=double> 
        struct MatlabVS { 
            // Prevent constructors 
            NO_CONSTRUCTORS(MatlabVS)

            // Setup the vector 
            typedef Optizelle::Matlab::Vector Vector; 

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
            // 0} where y > 0
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
            public Optizelle::ScalarValuedFunction <double,MatlabVS>
        {
        private:
            // Underlying function 
            std::unique_ptr <mxArrayPtr> data;

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(ScalarValuedFunction)

            // Create a function 
            template <typename A>
            ScalarValuedFunction(A && data_) :
                data(new A(std::move(data_)))
            {}

            // <- f(x) 
            double eval(Vector const & x) const; 

            // g = grad f(x) 
            void grad(Vector const & x,Vector & grad) const; 

            // H_dx = hess f(x) dx 
            void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)const;
        };

        // A simple vector valued function interface, f : X -> Y
        struct VectorValuedFunction :
            public Optizelle::VectorValuedFunction<double,MatlabVS,MatlabVS>
        {
        private:
            // Create some type shortcuts
            typedef Vector X_Vector; 
            typedef Vector Y_Vector; 

            // Name of this function
            std::string const name;

            // Underlying function 
            std::unique_ptr <mxArrayPtr> data;

        public:
            // Prevent constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(VectorValuedFunction)

            // Create a function 
            template <typename A>
            VectorValuedFunction(
                std::string const & name_,
                A && data_
            ) :
                name(name_), data(new A(std::move(data_)))
            {}

            // y=f(x)
            void eval(X_Vector const & x,Y_Vector & y) const;

            // y=f'(x)dx 
            void p(X_Vector const & x,X_Vector const & dx,Y_Vector & y) const;

            // xhat=f'(x)*dy
            void ps(X_Vector const & x,const Y_Vector &dy,X_Vector &xhat) const; 
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
            public Optizelle::Operator <double,MatlabVS,MatlabVS>
        {
        private:
            // Create some type shortcuts
            typedef Vector X_Vector; 
            typedef Vector Y_Vector; 

            // Name of this function
            std::string const name;

            // Underlying operator
            std::unique_ptr <mxArrayPtr> data;

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
            State <ProblemClass> & mxstate;
            typename ProblemClass::State::t const & state;

        public:
            // Prevent constructors
            NO_DEFAULT_COPY_ASSIGNMENT(Operator)

            // Create an operator
            template <typename A>
            Operator(
                std::string const & name_,
                A && data_,
                State <ProblemClass> & mxstate_,
                typename ProblemClass::State::t const & state_
            ) :
                name(name_),
                data(new A(std::move(data_))),
                mxstate(mxstate_),
                state(state_)
            {}

            // y = A(x)
            void eval(X_Vector const & x,Y_Vector & y) const {
                // Convert the state to a Matlab state
                mxstate.toMatlab(state);

                // Apply the operator to the state, x, and y
                auto eval = mxGetField(data->get(),0,"eval");
                auto ret_err = mxArray_CallObject2(
                    eval,
                    mxstate.data->get(),
                    x.data->get());
                
                // Check errors
                if(std::get <1> (ret_err))
                    throw Matlab::Exception::t( __LOC__
                        + ", evaluation of the eval function in the operator "
                        + name + " failed");
            
                // Copy the returned value into y 
                y.data.reset(new mxManaged(std::move(std::get<0> (ret_err))));
            }
        };
        
        // Converts elements from C++ to Matlab 
        namespace toMatlab {
            // Sets a real in a Matlab state 
            void Real(
                std::string const & name,
                double const & value,
                mxArrayPtr & mxstate
            );

            // Sets a natural in a Matlab state 
            void Natural(
                std::string const & name,
                Optizelle::Natural const & value,
                mxArrayPtr & mxstate
            );
           
            // Sets a parameter in a Matlab state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<mxManaged (enum_t const &)> const & toMatlab,
                enum_t const & value,
                mxArrayPtr & mxstate 
            ) {
                auto item = toMatlab(value);
                mxSetField(mxstate.get(),0,name.c_str(),item.release());
            }
            
            // Sets a vector in a Matlab state 
            void Vector(
                std::string const & name,
                Matlab::Vector const & value,
                mxArrayPtr & mxstate 
            );
            
            // Sets a list of vectors in a Matlab state 
            void VectorList(
                std::string const & name,
                std::list <Matlab::Vector> const & values,
                mxArrayPtr & mxstate 
            );
        
            // Sets restart vectors in Matlab 
            void Vectors(
                Matlab::Vectors const & values,
                mxArrayPtr & mxvalues 
            );

            // Sets restart reals in Matlab 
            void Reals(
                Matlab::Reals const & values,
                mxArrayPtr & mxvalues 
            );
            
            // Sets restart naturals in Matlab 
            void Naturals(
                Matlab::Naturals const & values,
                mxArrayPtr & mxvalues 
            );
            
            // Sets restart parameters in Matlab 
            void Params(
                Matlab::Params const & values,
                mxArrayPtr & mxvalues 
            );
        }

        // Converts elements from Matlab to C++ 
        namespace fromMatlab {
            // Sets a real in a C++ state 
            void Real(
                std::string const & name,
                mxArrayPtr const & mxstate,
                double & value
            );
            
            // Sets a natural in a C++ state 
            void Natural(
                std::string const & name,
                mxArrayPtr const & mxstate,
                Optizelle::Natural & value
            );
           
            // Sets a param C++ state 
            template <typename enum_t>
            void Param(
                std::string const & name,
                std::function<enum_t(mxArrayPtr const &)> const & fromMatlab,
                mxArrayPtr const & obj,
                enum_t & value
            ) {
                auto item = mxUnmanaged(mxGetField(obj.get(),0,name.c_str()));
                value = fromMatlab(item);
            }
            
            // Sets a vector in a C++ state 
            void Vector(
                std::string const & name,
                mxArrayPtr const & obj,
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
            template <typename ProblemClass>
            void ScalarValuedFunction(
                std::string const & name,
                Functions <ProblemClass> const & fns,
                std::unique_ptr <MxScalarValuedFunction> & value
            ) {
                value.reset(new Matlab::ScalarValuedFunction(
                    mxUnmanaged(mxGetField(fns.data->get(),0,name.c_str()))));
            }
            
            // Sets a vector-valued function in a C++ function bundle
            template <typename ProblemClass>
            void VectorValuedFunction(
                std::string const & name,
                Functions <ProblemClass> const & fns,
                std::unique_ptr <MxVectorValuedFunction> & value
            ) {
                value.reset(new Matlab::VectorValuedFunction(
                    name,
                    mxUnmanaged(mxGetField(fns.data->get(),0,name.c_str()))));
            }
            
            // Sets an operator in a C++ function bundle 
            template <typename ProblemClass>
            void Operator(
                std::string const & name,
                Functions <ProblemClass> const & mxfns,
                State <ProblemClass> & mxstate,
                typename ProblemClass::State::t const & state,
                std::unique_ptr <MxOperator> & value
            ) {
                value.reset(new Matlab::Operator <ProblemClass> (
                    name,
                    mxUnmanaged(mxGetField(mxfns.data->get(),0,name.c_str())),
                    mxstate,
                    state));
            }
        
            // Sets restart vectors in C++ 
            void Vectors(
                Matlab::Vector const & vec,
                mxArrayPtr const & mxvalues,
                Matlab::Vectors & values
            );
            
            // Sets restart reals in C++ 
            void Reals(
                mxArrayPtr const & mxvalues,
                Matlab::Reals & values
            );
            
            // Sets restart naturals in C++ 
            void Naturals(
                mxArrayPtr const & mxvalues,
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
                mxManaged mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxUnconstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                void toMatlab(
                    typename MxUnconstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxUnconstrained::State::t & state
                );
                void fromMatlab(
                    mxArrayPtr const & mxstate,
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
                    Matlab::Functions <ProblemClass> const & mxfns,
                    Matlab::State <ProblemClass> & mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxUnconstrained::Functions::t & fns 
                ) {
                    fromMatlab::ScalarValuedFunction("f",mxfns,fns.f);
                    fromMatlab::Operator <ProblemClass> (
                        "PH",mxfns,mxstate,state,fns.PH);
                }
                void fromMatlab(
                    Matlab::Functions <MxUnconstrained> const & mxfns,
                    Matlab::State <MxUnconstrained> & mxstate,
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
                mxManaged mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                void toMatlab(
                    typename MxEqualityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxEqualityConstrained::State::t & state
                );
                void fromMatlab(
                    mxArrayPtr const & mxstate,
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
                    Matlab::Functions <ProblemClass> const & mxfns,
                    Matlab::State <ProblemClass> & mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxEqualityConstrained::Functions::t & fns 
                ) {
                    fromMatlab::VectorValuedFunction("g",mxfns,fns.g);
                    fromMatlab::Operator <ProblemClass> ("PSchur_left",
                        mxfns,mxstate,state,fns.PSchur_left);
                    fromMatlab::Operator <ProblemClass> ("PSchur_right",
                        mxfns,mxstate,state,fns.PSchur_right);
                }
                void fromMatlab(
                    Matlab::Functions <MxEqualityConstrained> const & mxfns,
                    Matlab::State <MxEqualityConstrained> & mxstate,
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
                mxManaged mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab_(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                void toMatlab(
                    typename MxInequalityConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab_(
                    mxArrayPtr const & mxstate,
                    typename MxInequalityConstrained::State::t & state
                );
                void fromMatlab(
                    mxArrayPtr const & mxstate,
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
                    Matlab::Functions <ProblemClass> const & mxfns,
                    Matlab::State <ProblemClass> & mxstate,
                    typename ProblemClass::State::t const & state,
                    typename MxInequalityConstrained::Functions::t & fns 
                ) {
                    fromMatlab::VectorValuedFunction("h",mxfns,fns.h);
                }
                void fromMatlab(
                    Matlab::Functions <MxInequalityConstrained> const & mxfns,
                    Matlab::State <MxInequalityConstrained> & mxstate,
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
                mxManaged mxCreate();

                // Convert a C++ state to a Matlab state 
                void toMatlab(
                    typename MxConstrained::State::t const & state,
                    mxArrayPtr & mxstate
                );
                
                // Convert a Matlab state to C++ 
                void fromMatlab(
                    mxArrayPtr const & mxstate,
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
                    Matlab::Functions <MxConstrained> const & mxfns,
                    Matlab::State <MxConstrained> & mxstate,
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
