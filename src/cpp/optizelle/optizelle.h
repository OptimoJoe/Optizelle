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

#ifndef OPTIZELLE_H
#define OPTIZELLE_H

#include<list>
#include<vector>
#include<limits>
#include<cmath>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<memory>
#include<functional>
#include<algorithm>
#include<numeric>
#include<optizelle/linalg.h>

// Putting this into a class prevents its construction.  Essentially, we use
// this trick in order to create modules like in ML.  It also allows us to
// created templated namespaces.
#define NO_CONSTRUCTORS(Name) \
    Name() = delete; \
    Name(Name const &) = delete; \
    Name & operator = (Name const &) = delete; \
    ~Name() = delete;

// Disallows copying or assigning.  This is useful for things like the
// StateManipulator or FunctionModification classes.
#define NO_COPY_ASSIGNMENT(Name) \
    Name(Name const &) = delete; \
    Name & operator = (Name const &) = delete;

// Disallows copying, assigning, or default construction.  This is useful for
// for things like the State::t types.
#define NO_DEFAULT_COPY_ASSIGNMENT(Name) \
    Name() = delete; \
    Name(Name const &) = delete; \
    Name& operator = (Name const &) = delete;

namespace Optizelle{

    // A scalar valued function interface, f : X -> R
    template <
        typename Real,
        template <typename> class XX
    >
    struct ScalarValuedFunction {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

        // <- f(x) 
        virtual Real operator () (Vector const & x) const = 0;

        // grad = grad f(x) 
        virtual void grad(Vector const & x,Vector & grad) const = 0;

        // H_dx = hess f(x) dx 
        virtual void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)
            const = 0;

        // Allow a derived class to deallocate memory
        virtual ~ScalarValuedFunction() {}
    };
    
    template <
        typename Real,
        template <typename> class XX
    >
    struct ScalarValuedFunctionModifications {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

        // Disallow constructors
        NO_COPY_ASSIGNMENT(ScalarValuedFunctionModifications); 

        // Use an empty default constructor
        ScalarValuedFunctionModifications() {}

        // Allow derived classes to deallocate memory
        virtual ~ScalarValuedFunctionModifications() {}

        // Merit function additions to the objective
        virtual Real merit(Vector const & x,Real const & f_x) const {
            return f_x;
        }

        // Stopping condition modification of the gradient
        virtual void grad_stop(
            Vector const & x,
            Vector const & grad,
            Vector& grad_stop
        ) const {
            X::copy(grad,grad_stop);
        }

        // Diagnostic modification of the gradient
        virtual void grad_diag(
            Vector const & x,
            Vector const & grad,
            Vector& grad_diag
        ) const {
            X::copy(grad,grad_diag);
        }

        // Modification of the gradient when finding a trial step
        virtual void grad_step(
            Vector const & x,
            Vector const & grad,
            Vector& grad_step
        ) const {
            X::copy(grad,grad_step);
        }

        // Modification of the gradient for a quasi-Newton method 
        virtual void grad_quasi(
            Vector const & x,
            Vector const & grad,
            Vector& grad_quasi
        ) const {
            X::copy(grad,grad_quasi);
        }

        // Modification of the gradient when solving for the equality multiplier
        virtual void grad_mult(
            Vector const & x,
            Vector const & grad,
            Vector& grad_mult
        ) const {
            X::copy(grad,grad_mult);
        }

        // Modification of the Hessian-vector product when finding a trial step
        virtual void hessvec_step(
            Vector const & x,
            Vector const & dx,
            Vector const & H_dx,
            Vector& Hdx_step 
        ) const {
            X::copy(H_dx,Hdx_step);
        }
    };


    // A vector valued function interface, f : X -> Y
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class YY 
    >
    struct VectorValuedFunction {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector; 
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector; 

        // y=f(x)
        virtual void operator () (X_Vector const & x,Y_Vector & y) const = 0;

         // y=f'(x)dx 
         virtual void p(
             X_Vector const & x,
             X_Vector const & dx,
             Y_Vector & y
         ) const = 0;

         // z=f'(x)*dy
         virtual void ps(
             X_Vector const & x,
             Y_Vector const & dy,
             X_Vector & z
         ) const= 0;
         
         // z=(f''(x)dx)*dy
         virtual void pps(
             X_Vector const & x,
             X_Vector const & dx,
             Y_Vector const & dy,
             X_Vector & z
         ) const = 0;
         
         // Allow a derived class to deallocate memory
         virtual ~VectorValuedFunction() {}
    };

    // Defines how we output messages to the user
    struct Messaging {
        // Prints a message
        virtual void print(std::string const & msg) const;

        // Prints an error
        virtual void error(std::string const & msg) const;

        // Allow a derived class to deallocate memory
        virtual ~Messaging();
    };

    // Which algorithm class do we use
    namespace AlgorithmClass{
        enum t : Natural{
            TrustRegion,            // Trust-Region algorithms
            LineSearch,             // Line-search algorithms
            UserDefined             // User provides the iterate 
        };

        // Converts the algorithm class to a string
        std::string to_string(t const & algorithm_class);
        
        // Converts a string to an algorithm class 
        t from_string(std::string const & algorithm_class);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }

    // Reasons why we stop the algorithm
    namespace StoppingCondition{
        enum t : Natural{
            NotConverged,            // Algorithm did not converge
            RelativeGradientSmall,   // Relative gradient was sufficiently small
            RelativeStepSmall,       // Relative change in the step is small
            MaxItersExceeded,        // Maximum number of iterations exceeded
            InteriorPointInstability,// Instability in the interior point method
            UserDefined              // Some user defined stopping condition 
        };

        // Converts the stopping condition to a string 
        std::string to_string(t const & opt_stop);

        // Converts a string to a stopping condition
        t from_string(std::string const & opt_stop);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const; 
        };
    }

    // Various operators for both Hessian approximations and preconditioners
    namespace Operators{
        enum t : Natural{
            Identity,          // Identity approximation
            ScaledIdentity,    // Scaled identity approximation
            BFGS,              // BFGS approximation
            InvBFGS,           // Inverse BFGS approximation
            SR1,               // SR1 approximation
            InvSR1,            // Inverse SR1 approximation
            UserDefined        // User defined operator (such as the true
                               // Hessian for Newton's method)
        };
        
        // Converts the operator type to a string 
        std::string to_string(t const & op);
        
        // Converts a string to a operator 
        t from_string(std::string const & op);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }

    // Different kinds of search directions 
    namespace LineSearchDirection{
        enum t : Natural{
            SteepestDescent,          // SteepestDescent 
            FletcherReeves,           // Fletcher-Reeves CG
            PolakRibiere,             // Polak-Ribiere CG
            HestenesStiefel,          // HestenesStiefel CG
            BFGS,                     // Limited-memory BFGS 
            NewtonCG                  // Newton-CG
        };
        
        // Converts the line-search direction to a string 
        std::string to_string(t const & dir);
        
        // Converts a string to a line-search direction 
        t from_string(std::string const & dir);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }

    // Different sorts of line searches
    namespace LineSearchKind{
        enum t : Natural{
            Brents,           // Brent's minimization
            GoldenSection,    // Golden-section search 
            BackTracking,     // BackTracking search 
            TwoPointA,        // Barzilai and Borwein's method A
            TwoPointB         // Barzilai and Borwein's method B
        };
            
        // Converts the line-search kind to a string 
        std::string to_string(t const & kind);
        
        // Converts a string to a line-search kind 
        t from_string(std::string const & kind);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const; 
        };

        // Determine whether or not the line-search checks the sufficient
        // decrease condition.
        bool is_sufficient_decrease(t const & kind);
    }
    
    // Different points in the optimization algorithm
    namespace OptimizationLocation{
        enum t : Natural{
            // Occurs at the start of the optimization function 
            BeginningOfOptimization,

            // Occurs before the initial function and gradient evaluation 
            BeforeInitialFuncAndGrad,

            // Occurs after the initial function and gradient evaluation 
            AfterInitialFuncAndGrad,
            
            // Occurs just before the main optimization loop 
            BeforeOptimizationLoop,

            // Occurs just before we take the optimization step x+dx
            BeforeSaveOld,

            // Occurs just before we take the optimization step x+dx
            BeforeStep,

            // Occurs before we calculate our new step.
            BeforeGetStep,
            
            // Occurs during a user defined get step calculation.
            GetStep,

            // Occurs after we take the optimization step x+dx, but before
            // we calculate the gradient based on this new step.  In addition,
            // after this point we set the objective value, f_x, to be
            // f_xpdx.
            AfterStepBeforeGradient,

            // Occurs just after the gradient computation with the new
            // trial step
            AfterGradient,

            // Occurs before we update our quasi-Newton information. 
            BeforeQuasi,

            // Occurs after we update our quasi-Newton information. 
            AfterQuasi,

            // This occurs last in the optimization loop.  At this point,
            // we have already incremented our optimization iteration and
            // checked our stopping condition
            EndOfOptimizationIteration,

            // This occurs prior to the computation of the line search
            BeforeLineSearch,

            // This occurs after a rejected trust-region step
            AfterRejectedTrustRegion,

            // This occurs after a rejected line-search step
            AfterRejectedLineSearch,

            // This occurs prior to checking the predicted versus actual
            // reduction in a trust-region method.
            BeforeActualVersusPredicted,

            // This occurs at the end of a Krylov iteration
            EndOfKrylovIteration,

            // This occurs at the end of all optimization 
            EndOfOptimization
        };
            
        // Converts the optimization location to a string 
        std::string to_string(t const & loc);
        
        // Converts a string to a line-search kind 
        t from_string(std::string const & loc);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }
    
    // Different problem classes
    namespace ProblemClass{
        enum t : Natural{
            Unconstrained,         // Unconstrained optimization 
            EqualityConstrained,   // Equality constrained optimization 
            InequalityConstrained, // Inequality constrained optimization 
            Constrained            // Fully constrained optimization 
        };

        // Converts the problem class to a string
        std::string to_string(t const & problem_class);

        // Converts a string to a problem class 
        t from_string(std::string const & problem_class);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }
    
    // Different truncated Krylov solvers 
    namespace KrylovSolverTruncated{
        enum t : Natural{
            ConjugateDirection,         // Conjugate direction 
            MINRES                      // MINRES 
        };

        // Converts the problem class to a string
        std::string to_string(t const & truncated_krylov);

        // Converts a string to a problem class 
        t from_string(std::string const & truncated_krylov);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    };
    
    
    // Different kinds of interior point methods
    namespace InteriorPointMethod{
        enum t : Natural{
            PrimalDual,          // Standard primal-dual interior point method 
            PrimalDualLinked,    // A primal dual IPM, but the primal and dual
                                 // variables are kept in lock step.
            LogBarrier           // Primal log-barrier method 
        };
        
        // Converts the interior point method to a string 
        std::string to_string(t const & ipm);
        
        // Converts a string to an interior point method 
        t from_string(std::string const & ipm);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }
    
    // Different schemes for adjusting the interior point centrality 
    namespace CentralityStrategy{
        enum t : Natural{
            Constant,           // We keep sigma fixed at each iteration.
            StairStep,          // If the relative improvement in the
                                // interior point parameter does not exceed
                                // that of the gradient, we do a constant
                                // reduction.  Otherwise, we hold sigma
                                // constant.
            PredictorCorrector  // On odd iterations, sigma=1, on even, sigma=0.
        };
        
        // Converts the centrality strategy to a string
        std::string to_string(t const & cstrat);
        
        // Converts a string to the cstrat
        t from_string(std::string const & cstrat);

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (std::string const & name) const;
        };
    }

    // A collection of miscellaneous diagnostics that help determine errors.
    namespace Diagnostics {
        // Returns the smallest positive non-Nan number between the two. 
        // If both are NaN, it will return NaN.
        template <typename Real>
        Real get_smallest(const Real x,const Real y) {
            return (x < y) || (y != y) ? x : y;
        }

        // Performs a 4-point finite difference directional derivative on
        // a scalar valued function f : X->R.  In other words, <- f'(x)dx.  We
        // accomplish this by doing a finite difference calculation on f.
        template <
            typename Real,
            template <typename> class XX
        >
        Real directionalDerivative(
            ScalarValuedFunction<Real,XX> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            Real const & epsilon
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx(X::init(x));

            // f(x+eps dx)
            X::copy(x,x_op_dx);
            X::axpy(epsilon,dx,x_op_dx);
            Real obj_xpes=f(x_op_dx);

            // f(x-eps dx)
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            Real obj_xmes=f(x_op_dx);

            // f(x+2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(2.*epsilon),dx,x_op_dx);
            Real obj_xp2es=f(x_op_dx);

            // f(x-2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.*epsilon),dx,x_op_dx);
            Real obj_xm2es=f(x_op_dx);

            // Calculate the directional derivative and return it
            Real dd=(obj_xm2es-Real(8.)*obj_xmes+Real(8.)*obj_xpes-obj_xp2es)
                /(Real(12.)*epsilon);
            return dd;
        }
        
        // Performs a 4-point finite difference directional derivative on
        // the gradient of a scalar valued function f : X->R.  In other words,
        // dd ~= hess f(x) dx.  We accomplish this by doing a finite difference
        // calculation on G where G(x)=grad f(x).
        template <
            typename Real,
            template <typename> class XX
        >
        void directionalDerivative(
            ScalarValuedFunction<Real,XX> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            Real const & epsilon,
            typename XX <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx(X::init(x));

            // Create an element to store the gradient at this point 
            X_Vector fgrad_x_op_dx(X::init(x));

            // Zero out the directional derivative
            X::zero(dd);

            // grad f(x+eps dx)
            X::copy(x,x_op_dx);
            X::axpy(epsilon,dx,x_op_dx);
            f.grad(x_op_dx,fgrad_x_op_dx);
            X::axpy(Real(8.),fgrad_x_op_dx,dd);

            // grad f(x-eps dx)
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            f.grad(x_op_dx,fgrad_x_op_dx);
            X::axpy(Real(-8.),fgrad_x_op_dx,dd);

            // grad f(x+2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(2.)*epsilon,dx,x_op_dx);
            f.grad(x_op_dx,fgrad_x_op_dx);
            X::axpy(Real(-1.),fgrad_x_op_dx,dd);

            // grad f(x-2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.)*epsilon,dx,x_op_dx);
            f.grad(x_op_dx,fgrad_x_op_dx);
            X::axpy(Real(1.),fgrad_x_op_dx,dd);

            // Finish the finite difference calculation 
            X::scal(Real(1.)/(Real(12.)*epsilon),dd);
        }

        // Performs a 4-point finite difference directional derivative on
        // a vector valued function f : X->Y. In other words, dd ~= f'(x)dx.
        // We accomplish this by doing a finite difference calculation on f.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        void directionalDerivative(
            VectorValuedFunction<Real,XX,YY> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            Real const & epsilon,
            typename YY <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx(X::init(x));

            // Create an element for f(x+eps dx), etc.
            Y_Vector f_x_op_dx(Y::init(dd));
            
            // Zero out the directional derivative
            Y::zero(dd);

            // f(x+eps dx)
            X::copy(x,x_op_dx);
            X::axpy(epsilon,dx,x_op_dx);
            f(x_op_dx,f_x_op_dx);
            Y::axpy(Real(8.),f_x_op_dx,dd);

            // f(x-eps dx)
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            f(x_op_dx,f_x_op_dx);
            Y::axpy(Real(-8.),f_x_op_dx,dd);

            // f(x+2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(2.)*epsilon,dx,x_op_dx);
            f(x_op_dx,f_x_op_dx);
            Y::axpy(Real(-1.),f_x_op_dx,dd);

            // f(x-2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.)*epsilon,dx,x_op_dx);
            f(x_op_dx,f_x_op_dx);
            Y::axpy(Real(1.),f_x_op_dx,dd);

            // Finish the finite difference calculation 
            Y::scal(Real(1.)/(Real(12.)*epsilon),dd);
        }
        
        // Performs a 4-point finite difference directional derivative on
        // the second derivative-adjoint of a vector valued function. In other
        // words, dd ~= (f''(x)dx)*dy.  In order to calculate this, we do a
        // finite difference approximation using g(x)=f'(x)*dy.  Therefore,
        // the error in the approximation should be in the dx piece.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        void directionalDerivative(
            VectorValuedFunction<Real,XX,YY> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            typename YY <Real>::Vector const & dy,
            Real const & epsilon,
            typename XX <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx(X::init(x));

            // Create an element for f'(x+eps dx)*dy, etc.
            X_Vector fps_xopdx_dy(X::init(dd));

            // Zero out the directional derivative
            X::zero(dd);

            // f'(x+eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            X::axpy(Real(8.),fps_xopdx_dy,dd);

            // f'(x-eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            X::axpy(Real(-8.),fps_xopdx_dy,dd);

            // f'(x+2 eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(Real(2.)*epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            X::axpy(Real(-1.),fps_xopdx_dy,dd);

            // f'(x-2 eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.)*epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            X::axpy(Real(1.),fps_xopdx_dy,dd);

            // Finish the finite difference calculation 
            X::scal(Real(1.)/(Real(12.)*epsilon),dd);
        }

        // Performs a finite difference test on the gradient of f where  
        // f : X->R is scalar valued.  In other words, we check grad f using f
        // and return the smallest relative error. 
        template <
            typename Real,
            template <typename> class XX
        >
        Real gradientCheck(
            Messaging const & msg,
            ScalarValuedFunction<Real,XX> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Calculate the gradient at the point x
            X_Vector f_grad(X::init(x));
            f.grad(x,f_grad);

            // Begin by calculating the directional derivative via the gradient
            Real dd_grad=X::innr(f_grad,dx);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the gradient.");
            Real min_rel_err(std::numeric_limits<Real>::quiet_NaN());
            for(Integer i=-2;i<=5;i++){
                Real epsilon=pow(Real(.1),int(i));
                Real dd=directionalDerivative <> (f,x,dx,epsilon);

                // Calculate the relative error
                Real rel_err=fabs(dd_grad-dd)
                    / (std::numeric_limits <Real>::epsilon()+fabs(dd_grad));

                // Calculate the smallest relative error seen so far 
                min_rel_err=get_smallest <> (rel_err,min_rel_err);

                // Print out the relative error
                std::stringstream ss;
                if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str());
            }
            
            // Return the function's smallest relative error
            return min_rel_err;
        }
        
        // Performs a finite difference test on the hessian of f where f : X->R
        // is scalar valued.  In other words, we check hess f dx using grad f.
        template <
            typename Real,
            template <typename> class XX
        >
        Real hessianCheck(
            Messaging const & msg,
            ScalarValuedFunction<Real,XX> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for the residual between the directional 
            // derivative computed Hessian-vector product and the true 
            // Hessian-vector product.
            X_Vector res(X::init(x));

            // Calculate hess f in the direction dx.  
            X_Vector hess_f_dx(X::init(x));
            f.hessvec(x,dx,hess_f_dx);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the Hessian.");
            Real min_rel_err(std::numeric_limits<Real>::quiet_NaN());
            for(Integer i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),int(i));
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),hess_f_dx,res);

                // Determine the relative error
                Real rel_err=sqrt(X::innr(res,res))
                    / (std::numeric_limits <Real>::epsilon()
                    + sqrt(X::innr(hess_f_dx,hess_f_dx)));

                // Calculate the smallest relative error seen so far 
                min_rel_err=get_smallest <> (rel_err,min_rel_err);

                // Print out the differences
                std::stringstream ss;
                if(i<0)ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str());
            }
            
            // Return the function's smallest relative error
            return min_rel_err;
        }
        
        // This tests the symmetry of the Hessian.  We accomplish this by
        // comparing <H(x)dx,dxx> to <dx,H(x)dxx>.
        template <
            typename Real,
            template <typename> class XX
        >
        Real hessianSymmetryCheck(
            Messaging const & msg,
            ScalarValuedFunction<Real,XX> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            typename XX <Real>::Vector const & dxx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Calculate hess f in the direction dx.  
            X_Vector H_x_dx(X::init(x));
            f.hessvec(x,dx,H_x_dx);
            
            // Calculate hess f in the direction dxx.  
            X_Vector H_x_dxx(X::init(x));
            f.hessvec(x,dxx,H_x_dxx);
            
            // Calculate <H(x)dx,dxx>
            Real innr_Hxdx_dxx = X::innr(H_x_dx,dxx);
            
            // Calculate <dx,H(x)dxx>
            Real innr_dx_Hxdxx = X::innr(dx,H_x_dxx);

            // Determine the absolute difference between the two.  This really
            // should be zero.
            Real diff=fabs(innr_Hxdx_dxx-innr_dx_Hxdxx);

            // Send a message with the result
            msg.print("Symmetry test on the Hessian of a scalar valued "
                "function.");
            std::stringstream ss;
            ss<< "The absolute err. between <H(x)dx,dxx> and <dx,H(x)dxx>: "
                << std::scientific << std::setprecision(16) << diff;
            msg.print(ss.str());
            
            // Return the absolute error in symmetry 
            return diff;
        }

        // Performs a finite difference test on the derivative of a
        // vector-valued function f.  Specifically, we check f'(x)dx using f.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        Real derivativeCheck(
            Messaging const & msg,
            VectorValuedFunction<Real,XX,YY> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            typename YY <Real>::Vector const & y
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            Y_Vector res(Y::init(y));

            // Calculate f'(x)dx 
            Y_Vector fp_x_dx(Y::init(y));
            f.p(x,dx,fp_x_dx);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the derivative of a "
                "vector-valued function.");
            Real min_rel_err(std::numeric_limits<Real>::quiet_NaN());
            for(Integer i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),int(i));
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                Y::axpy(Real(-1.),fp_x_dx,res);

                // Determine the relative error
                Real rel_err=sqrt(Y::innr(res,res))
                    / (std::numeric_limits <Real>::epsilon()
                    + sqrt(Y::innr(fp_x_dx,fp_x_dx)));
                
                // Calculate the smallest relative error seen so far 
                min_rel_err=get_smallest <> (rel_err,min_rel_err);

                // Print out the differences
                std::stringstream ss;
                if(i<0)ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str());
            }
            
            // Return the function's smallest relative error
            return min_rel_err; 
        }

        // Performs an adjoint check on the first-order derivative of a vector
        // valued function.  In other words, we check that
        // <f'(x)dx,dy> = <dx,f'(x)*dy>
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        Real derivativeAdjointCheck(
            Messaging const & msg,
            VectorValuedFunction<Real,XX,YY> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            typename YY <Real>::Vector const & dy
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Calculate f'(x)dx 
            Y_Vector fp_x_dx(Y::init(dy));
            f.p(x,dx,fp_x_dx);
            
            // Calculate f'(x)*dy 
            X_Vector fps_x_dy(X::init(dx));
            f.ps(x,dy,fps_x_dy);

            // Calculate <f'(x)dx,dy>
            Real innr_fpxdx_dy = Y::innr(fp_x_dx,dy);

            // Calculate <dx,f'(x)*dy>
            Real innr_dx_fpsxdy = X::innr(dx,fps_x_dy);

            // Determine the absolute difference between the two.  This really
            // should be zero.
            Real diff=fabs(innr_fpxdx_dy-innr_dx_fpsxdy);

            // Send a message with the result
            msg.print("Adjoint test on the first derivative of a vector "
                "valued function.");
            std::stringstream ss;
            ss<<"The absolute err. between <f'(x)dx,dy> and <dx,f'(x)*dy>: "
                << std::scientific << std::setprecision(16) << diff;
            msg.print(ss.str());
            
            // Return the absolute error in symmetry 
            return diff;
        }

        // Performs a finite difference test on the second-derivative-adjoint 
        // of a vector-valued function f.  Specifically, we check
        // (f''(x)dx)*dy using f'(x)*dy.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        Real secondDerivativeCheck(
            Messaging const & msg,
            VectorValuedFunction<Real,XX,YY> const & f,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            typename YY <Real>::Vector const & dy
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            X_Vector res(X::init(x));

            // Calculate (f''(x)dx)*dy
            X_Vector fpps_x_dx_dy(X::init(dx));
            f.pps(x,dx,dy,fpps_x_dx_dy);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the 2nd-derivative adj. "
                "of a vector-valued function.");
            Real min_rel_err(std::numeric_limits<Real>::quiet_NaN());
            for(Integer i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),int(i));
                directionalDerivative <> (f,x,dx,dy,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),fpps_x_dx_dy,res);

                // Determine the relative error
                Real rel_err=sqrt(X::innr(res,res))
                    / (std::numeric_limits <Real>::epsilon()
                    + sqrt(X::innr(fpps_x_dx_dy,fpps_x_dx_dy)));
                
                // Calculate the smallest relative error seen so far 
                min_rel_err=get_smallest <> (rel_err,min_rel_err);

                // Print out the differences
                std::stringstream ss;
                if(i<0)ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str());
            }
            
            // Return the function's smallest relative error
            return min_rel_err; 
        }
    }

    // A function that has free reign to manipulate or analyze the state.
    // This should be used cautiously.
    template <typename ProblemClass>
    struct StateManipulator {
        // Disallow constructors
        NO_COPY_ASSIGNMENT(StateManipulator);

        // Default constructor
        StateManipulator() {}
        
        // Application
        virtual void operator () (
            typename ProblemClass::Functions::t const & fns,
            typename ProblemClass::State::t & state,
            OptimizationLocation::t const & loc
        ) const {};

        // Allow the derived class to deallocate memory
        virtual ~StateManipulator() {}
    };
   
    // A state manipulator that's been customized in order to print diagonistic
    // information
    template <typename ProblemClass>
    struct DiagnosticManipulator : public StateManipulator <ProblemClass> {
    private:
        // A reference to an existing state manipulator 
        StateManipulator <ProblemClass> const & smanip;

        // A reference to the messsaging object
        Messaging const & msg;

    public:
        // Disallow constructors
        NO_COPY_ASSIGNMENT(DiagnosticManipulator);

        // Create a reference to an existing manipulator 
        explicit DiagnosticManipulator(
            StateManipulator <ProblemClass> const & smanip_,
            Messaging const & msg_
        ) : smanip(smanip_), msg(msg_) {}

        // Application
        void operator () (
            typename ProblemClass::Functions::t const & fns,
            typename ProblemClass::State:: t& state,
            OptimizationLocation::t const & loc
        ) const {

            // Call the internal manipulator 
            smanip(fns,state,loc);

            // Create some shortcuts
            Natural const & msg_level=state.msg_level;

            switch(loc){

            // Output the headers for the diagonstic information
            case OptimizationLocation::BeforeOptimizationLoop:
                if(msg_level >= 1) {
                    // Get the headers 
                    std::list <std::string> out;
                    ProblemClass::Printer::getStateHeader(state,out);
                    if(msg_level >= 2)
                        ProblemClass::Printer::getKrylovHeader(state,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                }
            // Output the overall state at the end of the optimization
            // iteration
            case OptimizationLocation::EndOfOptimizationIteration: 
            case OptimizationLocation::AfterRejectedTrustRegion:
            case OptimizationLocation::AfterRejectedLineSearch:
                if(msg_level >= 1) {
                    // Get the diagonstic information
                    std::list <std::string> out;

                    // Don't print out the iteration information if
                    // we reject the step
                    if(    loc==OptimizationLocation
                            ::AfterRejectedTrustRegion
                        || loc==OptimizationLocation
                            ::AfterRejectedLineSearch
                    )
                        ProblemClass::Printer
                            ::getState(fns,state,false,true,out);
                    else 
                        ProblemClass::Printer
                            ::getState(fns,state,false,false,out);

                    // Print out blank Krylov information 
                    ProblemClass::Printer::getKrylov(state,true,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                }
                break;

            // Output information at the end of each Krylov iteration
            case OptimizationLocation::EndOfKrylovIteration:
                if(msg_level >= 2) {
                    // Get the diagonstic information
                    std::list <std::string> out;

                    // Print out blank state information
                    ProblemClass::Printer::getState(fns,state,true,false,out);

                    // Print out the Krylov information 
                    ProblemClass::Printer::getKrylov(state,false,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                }
                break;

            default:
                break;
            }
        }
    };

    // This converts one manipulator to another.  In theory, the dynamic
    // casting can fail, so make sure to only use this when compatibility
    // can be guaranteed.
    template <typename Internal,typename External> 
    struct ConversionManipulator : public StateManipulator <External> {
    private:
        // A reference to the user-defined state manipulator
        StateManipulator <Internal> const & smanip;

    public:
        // Disallow constructors
        NO_COPY_ASSIGNMENT(ConversionManipulator);

        // Grab a copy of the internal manipulator
        explicit ConversionManipulator(
            const StateManipulator <Internal>& smanip_
        ) : smanip(smanip_) {}

        // Application
        void operator () (
            const typename External::Functions::t& fns_,
            typename External::State::t& state_,
            OptimizationLocation::t const & loc
        ) const {
            const typename Internal::Functions::t& fns
                =dynamic_cast <typename Internal::Functions::t const &> (fns_);
            typename Internal::State::t& state 
                =dynamic_cast <typename Internal::State::t &> (state_);
            smanip(fns,state,loc);
        }
    };

    // A series of utiilty functions used by the routines below.
    namespace Utility {
        // Checks whether all the labels in the list labels are actually
        // labels stored in the is_label function.  If not, this function
        // throws an error.
        void checkLabels(
            Messaging const & msg,
            std::function<bool(std::string const &)> const & is_label,
            std::list <std::string> const & labels, 
            std::string const & kind
        );

        // Combines two strings in a funny sort of way.  Basically, given
        // a and b, if a is empty, this function returns b.  However, if
        // a is nonempty, it returns a.
        struct combineStrings : public std::binary_function
            <std::string const &,std::string const &,std::basic_string <char> >
        {
            std::basic_string <char> operator() (
                std::string const & a,std::string const & b) const;
        };
      
        // This checks the parameters and prints and error in there's a problem.
        void checkParams(
            Messaging const & msg,
            std::function
                <std::string(std::string const &,std::string const &)> const &
                checkParamVal,
            std::pair <std::list <std::string>,std::list <std::string> > const &
                params 
        );

        // Converts a variety of basic datatypes to strings
        std::ostream& formatReal(std::ostream& out);
        std::ostream& formatInt(std::ostream& out);
        std::ostream& formatString(std::ostream& out); 

        // Converts anything to a string.
        std::string atos(const double& x);
        std::string atos(Natural const & x);
        std::string atos(std::string const & x);
        std::string atos(KrylovStop::t const & x);

        // Blank separator for printing
        std::string const blankSeparator = ".         ";
    }
       
    // Routines that manipulate and support problems of the form
    // 
    // min_{x \in X} f(x)
    //
    // where f : X -> R.
    template <
        typename Real,
        template <typename> class XX
    > 
    struct Unconstrained {
    protected:
        // Different nonlinear-CG directions 
        struct NonlinearCGDirections {
            enum t : Natural{
                HestenesStiefel,        
                PolakRibiere,           
                FletcherReeves          
            };
        };
   
        // Return whether the line-search hit the bounds when it
        // terminated.  This is important to determine if we need
        // to find a different area to bracket.
        struct LineSearchTermination{
            enum t : Natural{
                Min,
                Max,
                Between
            };
        };

    public:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        
        // Disallow constructors
        NO_CONSTRUCTORS(Unconstrained);

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State);
                
            // Internal state of the optimization
            struct t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t);

                // ------------- GENERIC ------------- 

                // Tolerance for the gradient stopping condition
                Real eps_grad;

                // Tolerance for the step length stopping criteria
                Real eps_dx;

                // Number of control objects to store in a quasi-Newton method
                Natural stored_history;

                // Number of failed iterations before we reset the history for
                // quasi-Newton methods
                Natural history_reset;

                // Current iteration
                Natural iter;

                // Maximum number of optimization iterations
                Natural iter_max;

                // Why we've stopped the optimization
                StoppingCondition::t opt_stop;

                // Current number of Krylov iterations taken
                Natural krylov_iter;

                // Maximum number of iterations in the Krylov method
                Natural krylov_iter_max;

                // Total number of Krylov iterations taken
                Natural krylov_iter_total;

                // The maximum number of vectors we orthogonalize against in 
                // the Krylov method.  For something like CG, this is 1.
                Natural krylov_orthog_max;

                // Why the Krylov method was last stopped
                KrylovStop::t krylov_stop;

                // Relative error in the Krylov method
                Real krylov_rel_err;

                // Stopping tolerance for the Krylov method
                Real eps_krylov;

                // Truncated Krylov solver
                KrylovSolverTruncated::t krylov_solver;

                // Algorithm class
                AlgorithmClass::t algorithm_class;

                // Preconditioner for the Hessian
                Operators::t PH_type;

                // Hessian approximation
                Operators::t H_type;

                // Norm of a typical tradient
                Real norm_gradtyp;

                // Norm of a typical trial step
                Real norm_dxtyp;

                // Optimization variable 
                X_Vector x; 
                
                // Gradient, possibly of the objective, possibly of the
                // Lagrangian.  It depends on the context.
                X_Vector grad;
                
                // Trial step 
                X_Vector dx;
                
                // Old optimization variable 
                X_Vector x_old; 
                
                // Old gradient 
                X_Vector grad_old;
                
                // Old trial step 
                X_Vector dx_old;

                // Difference in prior gradients
                std::list <X_Vector> oldY;

                // Difference in prior steps
                std::list <X_Vector> oldS;

                // Current value of the objective function 
                Real f_x;

                // Objective function at the trial step
                Real f_xpdx;

                // Messaging level
                Natural msg_level;
                
                // ------------- TRUST-REGION ------------- 

                // Trust region radius
                Real delta;

                // Trust-region parameter for checking whether a step has been
                // accepted
                Real eta1;

                // Trust-region parameter for checking whether we enlarge 
                // the trust-region radius 
                Real eta2;

                // Actual reduction 
                Real ared;

                // Predicted reduction
                Real pred;

                // Number of rejected trust-region steps
                Natural rejected_trustregion;

                // ------------- LINE-SEARCH ------------- 

                // Base line-search step length
                Real alpha0;
                
                // Actual line-search step length
                Real alpha;

                // Parameter that helps govern the sufficient decrease
                Real c1;

                // Current number of iterations used in the line-search
                Natural linesearch_iter;

                // Maximum number of iterations used in the line-search
                Natural linesearch_iter_max;

                // Total number of line-search iterations computed
                Natural linesearch_iter_total;

                // Stopping tolerance for the line-search
                Real eps_ls;

                // Search direction type
                LineSearchDirection::t dir;

                // Type of line-search 
                LineSearchKind::t kind;

                // Initialization constructors
                explicit t(X_Vector const & x_) :
                    eps_grad(1e-8),
                    eps_dx(1e-8),
                    stored_history(0),
                    history_reset(5),
                    iter(1),
                    iter_max(10),
                    opt_stop(StoppingCondition::NotConverged),
                    krylov_iter(0),
                    krylov_iter_max(10),
                    krylov_iter_total(0),
                    krylov_orthog_max(1),
                    krylov_stop(KrylovStop::RelativeErrorSmall),
                    krylov_rel_err(0.),
                    eps_krylov(1e-2),
                    krylov_solver(KrylovSolverTruncated::ConjugateDirection),
                    algorithm_class(AlgorithmClass::TrustRegion),
                    PH_type(Operators::Identity),
                    H_type(Operators::UserDefined),
                    norm_gradtyp(std::numeric_limits<Real>::quiet_NaN()),
                    norm_dxtyp(std::numeric_limits<Real>::quiet_NaN()),
                    x(X::init(x_)),
                    grad(X::init(x_)),
                    dx(X::init(x_)),
                    x_old(X::init(x_)),
                    grad_old(X::init(x_)),
                    dx_old(X::init(x_)),
                    oldY(),
                    oldS(), 
                    f_x(std::numeric_limits<Real>::quiet_NaN()),
                    f_xpdx(std::numeric_limits<Real>::quiet_NaN()),
                    msg_level(1),
                    delta(1.),
                    eta1(.1),
                    eta2(.9),
                    ared(1.),
                    pred(1.),
                    rejected_trustregion(0),
                    alpha0(1.),
                    alpha(1.),
                    c1(1e-4),
                    linesearch_iter(0),
                    linesearch_iter_max(5),
                    linesearch_iter_total(0),
                    eps_ls(1e-2),
                    dir(LineSearchDirection::SteepestDescent),
                    kind(LineSearchKind::GoldenSection)
                {
                    X::copy(x_,x);
                }
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // Check that we have a valid set of parameters.  
            static void check_(Messaging const & msg,t const & state) {
                   
                // Use this to build an error message
                std::stringstream ss;
                
                // Check that the tolerance for the gradient stopping condition
                // is positive
                if(state.eps_grad <= Real(0.)) 
                    ss << "The tolerance for the gradient stopping condition "
                        "must be positive: eps_grad = " << state.eps_grad;
            
                // Check that the tolerance for the step length stopping
                // condition is positive
                else if(state.eps_dx <= Real(0.)) 
                    ss << "The tolerance for the step length stopping "
                        "condition must be positive: eps_dx = " << state.eps_dx;
        
                // Check that the current iteration is positive
                else if(state.iter == 0) 
                    ss << "The current optimization iteration must be "
                        "positive: iter = " << state.iter;

                // Check that the maximum iteration is positive
                else if(state.iter_max == 0) 
                    ss << "The maximum optimization iteration must be "
                        "positive: iter_max = " << state.iter_max;

                // Check that the maximum Krylov iteration is positive
                else if(state.krylov_iter_max == 0) 
                    ss << "The maximum Krylov iteration must be "
                        "positive: krylov_iter_max = " << state.krylov_iter_max;

                // Check that the number of vectors we orthogonalize against
                // is at least 1.
                else if(state.krylov_orthog_max == 0) 
                    ss << "The maximum number of vectors the Krylov method"
                    "orthogonalizes against must be positive: "
                    "krylov_orthog_max = " << state.krylov_orthog_max;

                // Check that relative error in the Krylov method is nonnegative
                else if(state.krylov_rel_err < Real(0.)) 
                    ss << "The relative error in the Krylov method must be "
                        "nonnegative: krylov_rel_err = " <<state.krylov_rel_err;
                
                // Check that the stopping tolerance for the Krylov method is
                // positive
                else if(state.eps_krylov <= Real(0.)) 
                    ss << "The tolerance for the Krylov method stopping "
                        "condition must be positive: eps_krylov = "
                    << state.eps_krylov;

                // Check that the norm of a typical gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_gradtyp < Real(0.)
                    || (state.iter!=1 && state.norm_gradtyp!=state.norm_gradtyp)
                ) 
                    ss << "The norm of a typical gradient must be nonnegative: "
                        "norm_gradtyp = " << state.norm_gradtyp; 

                // Check that the norm of a typical trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_dxtyp < Real(0.)
                    || (state.iter!=1 && state.norm_dxtyp!=state.norm_dxtyp)
                ) 
                    ss << "The norm of a typical trial step must be "
                        "nonnegative: norm_dxtyp = " << state.norm_dxtyp; 

                // Check that the objective value isn't a NaN past
                // iteration 1
                else if(state.iter!=1 && state.f_x!=state.f_x)
                    ss<< "The objective value must be a number: f_x = "
                        << state.f_x;

                // Check that the objective value at a trial step isn't
                // a NaN past iteration 1
                else if(state.iter!=1
                     && state.f_xpdx != state.f_xpdx
                ) 
                    ss << "The objective value at the trial step must be a "
                        "number: f_xpdx = " << state.f_xpdx;

                // Check that the trust-region radius is nonnegative 
                else if(state.delta<Real(0.))
                    ss<< "The trust-region radius must be nonnegative: delta = "
                        << state.delta; 

                // Check that the predicted vs. actual reduction tolerance
                // is between 0 and 1
                else if(state.eta1 < Real(0.) || state.eta1 > Real(1.))
                    ss << "The tolerance for whether or not we accept a "
                        "trust-region step must be between 0 and 1: eta1 = "
                        << state.eta1;
                
                // Check that the other predicted vs. actual reduction tolerance
                // is between 0 and 1
                else if(state.eta2 < Real(0.) || state.eta2 > Real(1.))
                    ss << "The tolerance for whether or not we increase the "
                        "trust-region radius must be between 0 and 1: eta2 = "
                        << state.eta2;

                // Check that eta2 > eta1
                else if(state.eta1 >= state.eta2) 
                    ss << "The trust-region tolerances for accepting steps "
                        "must satisfy the relationship that eta1 < eta2: "
                        "eta1 = " << state.eta1 << ", eta2 = " << state.eta2;

                // Check that the base line-search step length is positive 
                else if(state.alpha0 <= Real(0.)) 
                    ss << "The base line-search step length must be positive: "
                        "alpha0 = " << state.alpha0;

                // Check that the sufficient decrease parameter lies between
                // 0 and 1.
                else if(state.c1 <= Real(0.) || state.c1 >= Real(1.)) 
                    ss << "The sufficient decrease parameter must lie between "
                        "0 and 1: c1 = " << state.c1;

                // Check that the stopping tolerance for the line-search
                // methods is positive
                else if(state.eps_ls <= Real(0.)) 
                    ss << "The tolerance for the line-search stopping "
                        "condition must be positive: eps_ls = " << state.eps_ls;

                // Check that the number of line-search iterations is positve 
                else if(state.linesearch_iter_max <= 0) 
                    ss << "The maximum number of line-search iterations must "
                        "be positive: linesearch_iter_max = "
                        << state.linesearch_iter_max;

                // If we're using a golden-section search, make sure we allow
                // at least two iterations.
                else if(state.kind==LineSearchKind::GoldenSection && 
                    state.linesearch_iter_max <= 1
                ) 
                    ss << "When using a golden-seciton search, we require at "
                        "least 2 line-search iterations: linesearch_iter_max = "
                        << state.linesearch_iter_max;

                // If we're using the Barzilai-Borwein two point Hessian
                // approximation, make sure the direction is steepest descent.
                else if((state.kind==LineSearchKind::TwoPointA ||
                    state.kind==LineSearchKind::TwoPointB) &&
                    state.dir!=LineSearchDirection::SteepestDescent
                )
                    ss << "When using the Barzilai-Borwein two point Hessian "
                        "approximation line-search, the search direction must "
                        "be set to SteepestDescent: dir = "
                        << LineSearchDirection::to_string(state.dir);

                // If there's an error, print it
                if(ss.str()!="") msg.error(ss.str());
            }
            static void check(Messaging const & msg,t const & state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
            }
        };

        // Utilities for restarting the optimization
        struct Restart {
            // Disallow constructors
            NO_CONSTRUCTORS(Restart);

            // Holds restart information
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <Natural> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;

            // Checks whether we have a valid real label.
            struct is_real :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( name == "eps_grad" || 
                        name == "eps_dx" || 
                        name == "krylov_rel_err" || 
                        name == "eps_krylov" || 
                        name == "norm_gradtyp" || 
                        name == "norm_dxtyp" || 
                        name == "f_x" || 
                        name == "f_xpdx" ||
                        name == "delta" || 
                        name == "eta1" || 
                        name == "eta2" || 
                        name == "ared" || 
                        name == "pred" || 
                        name == "alpha0" || 
                        name == "alpha" || 
                        name == "c1" || 
                        name == "eps_ls"
                    ) 
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have a valid natural number label.
            struct is_nat :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( name == "stored_history" ||
                        name == "history_reset" || 
                        name == "iter" || 
                        name == "iter_max" || 
                        name == "krylov_iter" || 
                        name == "krylov_iter_max" ||
                        name == "krylov_iter_total" || 
                        name == "krylov_orthog_max" ||
                        name == "msg_level" ||
                        name == "rejected_trustregion" || 
                        name == "linesearch_iter" || 
                        name == "linesearch_iter_max" ||
                        name == "linesearch_iter_total" 
                    ) 
                        return true;
                    else
                        return false;
                }
            };
           
            // Checks whether we have a valid parameter label.
            struct is_param:
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( name == "krylov_solver" ||
                        name == "algorithm_class" || 
                        name == "opt_stop" || 
                        name == "krylov_stop" ||
                        name == "H_type" || 
                        name == "PH_type" ||
                        name == "dir" || 
                        name == "kind" 
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid variable label
            struct is_x : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( name == "x" || 
                        name == "grad" || 
                        name == "dx" || 
                        name == "x_old" || 
                        name == "grad_old" || 
                        name == "dx_old" || 
                        name.substr(0,5)=="oldY_" || 
                        name.substr(0,5)=="oldS_" 
                    ) 
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have valid labels
            static void checkLabels(
                Messaging const & msg,
                Reals const & reals,
                Nats const & nats,
                Params const & params,
                X_Vectors const & xs
            ) {
                Utility::checkLabels(msg,is_real(),reals.first," real name: ");
                Utility::checkLabels(msg,is_nat(),nats.first," natural name: ");
                Utility::checkLabels
                    (msg,is_param(),params.first," paramater name: ");
                Utility::checkLabels(msg,is_x(),xs.first," variable name: ");
            }

            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string const &,std::string const &,std::string>
            {
                std::string operator() (
                    std::string const & label,
                    std::string const & val
                ) const {
                    // Create a base message
                    const std::string base
                        ="During serialization, found an invalid ";

                    // Used to build the message 
                    std::stringstream ss;

                    // Check the truncated Krylov solver 
                    if(label=="krylov_solver") {
                        if(!KrylovSolverTruncated::is_valid()(val))
                            ss << base << "truncated Krylov solver: " << val;

                    // Check the algorithm class
                    } else if(label=="algorithm_class") {
                        if(!AlgorithmClass::is_valid()(val))
                            ss << base << "algorithm class: " << val;

                    // Check the optimization stopping conditions
                    } else if(label=="opt_stop"){
                        if(!StoppingCondition::is_valid()(val))
                            ss << base << "stopping condition: " << val;

                    // Check the Krylov stopping conditions
                    } else if(label=="krylov_stop"){
                        if(!KrylovStop::is_valid()(val)) 
                            ss <<base <<"Krylov stopping condition: " << val;

                    // Check the Hessian type
                    } else if(label=="H_type"){
                        if(!Operators::is_valid()(val))
                            ss << base<<"Hessian type: " << val;

                    // Check the type of the preconditioner
                    } else if(label=="PH_type"){
                        if(!Operators::is_valid()(val))
                            ss << base <<"Hessian preconditioner type: " << val;

                    // Check the line-search direction
                    } else if(label=="dir"){
                        if(!LineSearchDirection::is_valid()(val)) 
                            ss << base << "line-search direction: " << val;

                    // Check the kind of line-search
                    } else if(label=="kind"){
                        if(!LineSearchKind::is_valid()(val)) 
                            ss << base << "line-search kind: " << val;
                    }
                    return ss.str();
                }
            };
            
            // Copy out all variables.
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs
            ) {
                // Move the memory of all variables into the list 
                xs.first.emplace_back("x");
                xs.second.emplace_back(std::move(state.x));
                xs.first.emplace_back("grad");
                xs.second.emplace_back(std::move(state.grad));
                xs.first.emplace_back("dx");
                xs.second.emplace_back(std::move(state.dx));
                xs.first.emplace_back("x_old");
                xs.second.emplace_back(std::move(state.x_old));
                xs.first.emplace_back("grad_old");
                xs.second.emplace_back(std::move(state.grad_old));
                xs.first.emplace_back("dx_old");
                xs.second.emplace_back(std::move(state.dx_old));

                // Write out the quasi-Newton information with sequential names
                {Natural i=1;
                for(typename std::list<X_Vector>::iterator y=state.oldY.begin();
                    y!=state.oldY.end();
                    y=state.oldY.begin()
                ){
                    std::stringstream ss;
                    ss << "oldY_" << i;
                    xs.first.emplace_back(ss.str());
                    xs.second.splice(xs.second.end(),state.oldY,y);
                }}

                // Write out the quasi-Newton information with sequential names
                {Natural i=1;
                for(typename std::list<X_Vector>::iterator s=state.oldS.begin();
                    s!=state.oldS.end();
                    s=state.oldS.begin()
                ){
                    std::stringstream ss;
                    ss << "oldS_" << i;
                    xs.first.emplace_back(ss.str());
                    xs.second.splice(xs.second.end(),state.oldS,s);
                }}
            }
            
            // Copy out all non-variables.  This includes reals, naturals,
            // and parameters
            static void stateToScalars(
                typename State::t & state, 
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                
                // Copy in all the real numbers 
                reals.first.emplace_back("eps_grad");
                reals.second.emplace_back(state.eps_grad);
                reals.first.emplace_back("eps_dx");
                reals.second.emplace_back(state.eps_dx);
                reals.first.emplace_back("krylov_rel_err");
                reals.second.emplace_back(state.krylov_rel_err);
                reals.first.emplace_back("eps_krylov");
                reals.second.emplace_back(state.eps_krylov);
                reals.first.emplace_back("norm_gradtyp");
                reals.second.emplace_back(state.norm_gradtyp);
                reals.first.emplace_back("norm_dxtyp");
                reals.second.emplace_back(state.norm_dxtyp);
                reals.first.emplace_back("f_x");
                reals.second.emplace_back(state.f_x);
                reals.first.emplace_back("f_xpdx");
                reals.second.emplace_back(state.f_xpdx);
                reals.first.emplace_back("delta");
                reals.second.emplace_back(state.delta);
                reals.first.emplace_back("eta1");
                reals.second.emplace_back(state.eta1);
                reals.first.emplace_back("eta2");
                reals.second.emplace_back(state.eta2);
                reals.first.emplace_back("ared");
                reals.second.emplace_back(state.ared);
                reals.first.emplace_back("pred");
                reals.second.emplace_back(state.pred);
                reals.first.emplace_back("alpha0");
                reals.second.emplace_back(state.alpha0);
                reals.first.emplace_back("alpha");
                reals.second.emplace_back(state.alpha);
                reals.first.emplace_back("c1");
                reals.second.emplace_back(state.c1);
                reals.first.emplace_back("eps_ls");
                reals.second.emplace_back(state.eps_ls);

                // Copy in all the natural numbers
                nats.first.emplace_back("stored_history");
                nats.second.emplace_back(state.stored_history);
                nats.first.emplace_back("history_reset");
                nats.second.emplace_back(state.history_reset);
                nats.first.emplace_back("iter");
                nats.second.emplace_back(state.iter);
                nats.first.emplace_back("iter_max");
                nats.second.emplace_back(state.iter_max);
                nats.first.emplace_back("krylov_iter");
                nats.second.emplace_back(state.krylov_iter);
                nats.first.emplace_back("krylov_iter_max");
                nats.second.emplace_back(state.krylov_iter_max);
                nats.first.emplace_back("krylov_iter_total");
                nats.second.emplace_back(state.krylov_iter_total);
                nats.first.emplace_back("krylov_orthog_max");
                nats.second.emplace_back(state.krylov_orthog_max);
                nats.first.emplace_back("msg_level");
                nats.second.emplace_back(state.msg_level);
                nats.first.emplace_back("rejected_trustregion");
                nats.second.emplace_back(state.rejected_trustregion);
                nats.first.emplace_back("linesearch_iter");
                nats.second.emplace_back(state.linesearch_iter);
                nats.first.emplace_back("linesearch_iter_max");
                nats.second.emplace_back(state.linesearch_iter_max);
                nats.first.emplace_back("linesearch_iter_total");
                nats.second.emplace_back(state.linesearch_iter_total);

                // Copy in all the parameters
                params.first.emplace_back("krylov_solver");
                params.second.emplace_back(
                    KrylovSolverTruncated::to_string(state.krylov_solver));
                params.first.emplace_back("algorithm_class");
                params.second.emplace_back(
                    AlgorithmClass::to_string(state.algorithm_class));
                params.first.emplace_back("opt_stop");
                params.second.emplace_back(
                    StoppingCondition::to_string(state.opt_stop));
                params.first.emplace_back("krylov_stop");
                params.second.emplace_back(
                    KrylovStop::to_string(state.krylov_stop));
                params.first.emplace_back("H_type");
                params.second.emplace_back(
                    Operators::to_string(state.H_type));
                params.first.emplace_back("PH_type");
                params.second.emplace_back(
                    Operators::to_string(state.PH_type));
                params.first.emplace_back("dir");
                params.second.emplace_back(
                    LineSearchDirection::to_string(state.dir));
                params.first.emplace_back("kind");
                params.second.emplace_back(
                    LineSearchKind::to_string(state.kind));
            }

            // Copy in all variables.  This assumes that the quasi-Newton
            // information is being read in order.
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs
            ) {
                typename std::list <X_Vector>::iterator x=xs.second.begin();
                for(typename std::list <std::string>::iterator name
                        =xs.first.begin();
                    name!=xs.first.end();
                ) {
                    // Make a copy of the current iterators.  We use these
                    // to remove elements
                    typename std::list <std::string>::iterator name0 = name;
                    typename std::list <X_Vector>::iterator x0 = x;

                    // Increment our primary iterators 
                    name++; x++;

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name0=="x") 
                        state.x = std::move(*x0);
                    else if(*name0=="grad")
                        state.grad = std::move(*x0);
                    else if(*name0=="dx")
                        state.dx = std::move(*x0);
                    else if(*name0=="x_old")
                        state.x_old = std::move(*x0);
                    else if(*name0=="grad_old")
                        state.grad_old = std::move(*x0);
                    else if(*name0=="dx_old")
                        state.dx_old = std::move(*x0);
                    else if(name0->substr(0,5)=="oldY_")
                        state.oldY.emplace_back(std::move(*x0));
                    else if(name0->substr(0,5)=="oldS_")
                        state.oldS.emplace_back(std::move(*x0));
                }
            }

            // Copy in all non-variables.  This includes reals, naturals,
            // and parameters
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy in any reals 
                typename std::list <Real>::iterator real=reals.second.begin();
                for(std::list <std::string>::iterator name=reals.first.begin();
                    name!=reals.first.end();
                    name++,real++
                ){
                    if(*name=="eps_grad") state.eps_grad=*real;
                    else if(*name=="eps_dx") state.eps_dx=*real;
                    else if(*name=="krylov_rel_err") state.krylov_rel_err=*real;
                    else if(*name=="eps_krylov") state.eps_krylov=*real;
                    else if(*name=="norm_gradtyp") state.norm_gradtyp=*real;
                    else if(*name=="norm_dxtyp") state.norm_dxtyp=*real;
                    else if(*name=="f_x") state.f_x=*real;
                    else if(*name=="f_xpdx") state.f_xpdx=*real;
                    else if(*name=="delta") state.delta=*real;
                    else if(*name=="eta1") state.eta1=*real;
                    else if(*name=="eta2") state.eta2=*real;
                    else if(*name=="ared") state.ared=*real;
                    else if(*name=="pred") state.pred=*real;
                    else if(*name=="alpha0") state.alpha0=*real;
                    else if(*name=="alpha") state.alpha=*real;
                    else if(*name=="c1") state.c1=*real;
                    else if(*name=="eps_ls") state.eps_ls=*real;
                }
            
                // Next, copy in any naturals
                typename std::list <Natural>::iterator nat=nats.second.begin();
                for(std::list <std::string>::iterator name=nats.first.begin();
                    name!=nats.first.end();
                    name++,nat++
                ){
                    if(*name=="stored_history") state.stored_history=*nat;
                    else if(*name=="history_reset") state.history_reset=*nat;
                    else if(*name=="iter") state.iter=*nat;
                    else if(*name=="iter_max") state.iter_max=*nat;
                    else if(*name=="krylov_iter") state.krylov_iter=*nat;
                    else if(*name=="krylov_iter_max")
                        state.krylov_iter_max=*nat;
                    else if(*name=="krylov_iter_total")
                        state.krylov_iter_total=*nat;
                    else if(*name=="krylov_orthog_max")
                        state.krylov_orthog_max=*nat;
                    else if(*name=="msg_level")
                        state.msg_level=*nat;
                    else if(*name=="rejected_trustregion")
                        state.rejected_trustregion=*nat;
                    else if(*name=="linesearch_iter")
                        state.linesearch_iter=*nat;
                    else if(*name=="linesearch_iter_max")
                        state.linesearch_iter_max=*nat;
                    else if(*name=="linesearch_iter_total")
                        state.linesearch_iter_total=*nat;
                }
                    
                // Next, copy in any parameters 
                std::list <std::string>::iterator param=params.second.begin();
                for(std::list <std::string>::iterator name=params.first.begin();
                    name!=params.first.end();
                    name++,param++
                ){
                    if(*name=="krylov_solver")
                        state.krylov_solver
                            =KrylovSolverTruncated::from_string(*param);
                    else if(*name=="algorithm_class")
                        state.algorithm_class
                            =AlgorithmClass::from_string(*param);
                    else if(*name=="opt_stop")
                        state.opt_stop=StoppingCondition::from_string(*param);
                    else if(*name=="krylov_stop")
                        state.krylov_stop=KrylovStop::from_string(*param);
                    else if(*name=="H_type")
                        state.H_type=Operators::from_string(*param);
                    else if(*name=="PH_type")
                        state.PH_type=Operators::from_string(*param);
                    else if(*name=="dir")
                        state.dir=LineSearchDirection::from_string(*param);
                    else if(*name=="kind")
                        state.kind=LineSearchKind::from_string(*param);
                }
            }
            
            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::Restart::stateToVectors(state,xs);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::stateToScalars(state,reals,nats,params);
            }
            
            // Capture data from structures controlled by the user.  Note,
            // we don't sort the oldY and oldS based on the prefix.  In fact,
            // we completely ignore this information.  Therefore, this routine
            // really depends on oldY and oldS to have their elements inserted
            // into vars in order.  In other words, oldY_1 must come before
            // oldY_2, etc.
            static void capture(
                Messaging const & msg,
                typename State::t & state,
                X_Vectors & xs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs);

                // Check the strings used to represent parameters
                Utility::checkParams(msg,checkParamVal(),params);

                // Copy in the variables 
                Unconstrained <Real,XX>::Restart::vectorsToState(state,xs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::scalarsToState(state,reals,nats,params);

                // Check that we have a valid state 
                State::check(msg,state);
            }
        };

        // All the functions required by an optimization algorithm.  Note, this
        // routine owns the memory for these operations.  
        struct Functions {
            // Disallow constructors
            NO_CONSTRUCTORS(Functions);

            // Actual storage of the functions required
            struct t{
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t);
                
                // Objective function
                std::unique_ptr <ScalarValuedFunction <Real,XX> > f;

                // Objective function modifications
                std::unique_ptr <ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Preconditioner for the Hessian of the objective
                std::unique_ptr <Operator <Real,XX,XX> > PH;

                // Symmetric, positive definite operator that changes the
                // scaling of the trust-region.
                std::unique_ptr <Operator <Real,XX,XX> > TRS;

                // Initialize all of the pointers to null
                t() : f(nullptr), PH(nullptr), TRS(nullptr) {}
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // The identity operator 
            struct Identity : public Operator <Real,XX,XX> {
                void operator () (X_Vector const & dx,X_Vector & result) const{
                    X::copy(dx,result);
                }
            };

            // The scaled identity Hessian approximation.  Specifically, use use
            // || grad || / (2 delta) I where delta is the current size of the
            // trust-region.  This forces us into the trust-region at each
            // iteration.
            class ScaledIdentity : public Operator <Real,XX,XX> {
            private:
                // Objective modifications
                ScalarValuedFunctionModifications <Real,XX> const & f_mod;

                // Current iterate
                X_Vector const & x;

                // Gradient of the objective 
                X_Vector const & grad;

                // Maximum size of the trust-region radius
                Real const & delta;
            public:
                ScaledIdentity(
                    typename Functions::t const & fns,
                    typename State::t const & state
                ) : f_mod(*(fns.f_mod)),
                    x(state.x),
                    grad(state.grad),
                    delta(state.delta)
                {};

                void operator () (X_Vector const & dx,X_Vector & result) const{
                    // Determine the norm of the gradient
                    X_Vector grad_step(X::init(grad));
                        f_mod.grad_step(x,grad,grad_step);
                    Real norm_grad=sqrt(X::innr(grad_step,grad_step));

                    // Copy in the direction and scale it
                    X::copy(dx,result);
                    X::scal(norm_grad/(Real(2.)*delta),result);
                }
            };

            // The BFGS Hessian approximation.  Note, the formula we normally
            // see for BFGS denotes the inverse Hessian approximation.  This is
            // not the inverse, but the true Hessian approximation. 
            class BFGS : public Operator <Real,XX,XX> {
            private:
                // Messaging device in case the quasi-Newton information is bad
                Messaging const & msg;

                // Stored quasi-Newton information
                std::list<X_Vector> const & oldY;
                std::list<X_Vector> const & oldS;
            public:
                BFGS(
                    Messaging const & msg_,
                    typename State::t const & state
                ) : msg(msg_), oldY(state.oldY), oldS(state.oldS) {};

                // Operator interface
                /* It's not entirely clear to me what the best implementation
                for this method really is.  In the following implementation, we
                require an additional k work elements where k is the number of
                stored gradient and position differences.  It's possible to
                reduce this to 1 or 2, but we need to compute redundant
                information.  It's also possible to implementation the compact
                representation, see "Representations of quasi-Newton matrices
                and their use in limited memory methods" from Byrd, Nocedal,
                and Schnabel.  The problem with that algorithm is that is
                requires machinery such as linear system solves that we don't
                current have.  It also works much better with matrices or
                multivectors of data and we don't require the user to provide
                these abstractions. */
                void operator () (X_Vector const & dx, X_Vector & result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the BFGS Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work;
                    for(int i=0;i<oldY.size();i++)
                        work.emplace_back(std::move(X::init(dx)));

                    // If we have no vectors in our history, we return the
                    // direction
                    X::copy(dx,result);
                    if(oldY.size() == 0) return;

                    // As a safety check, insure that the inner product
                    // between all the (s,y) pairs is positive
                    typename std::list <X_Vector>::const_iterator y0
                        = oldY.begin();
                    typename std::list <X_Vector>::const_iterator s0
                        = oldS.begin();
                    while(y0!=oldY.end()){
                        Real inner_y_s=X::innr(*y0++,*s0++);
                        if(inner_y_s <= Real(0.))
                            msg.error("Detected a (s,y) pair in BFGS that "
                                "possesed a nonpositive inner product");
                    }

                    // Othwerwise, we copy all of the trial step differences
                    // into the work space
                    typename std::list <X_Vector>::iterator Bisj_iter
                        =work.begin();
                    typename std::list <X_Vector>::const_iterator sk_iter
                        =oldS.begin();
                    while(Bisj_iter!=work.end())
                        X::copy((*sk_iter++),(*Bisj_iter++));

                    // Keep track of the element Bisi
                    typename std::list <X_Vector>::const_iterator Bisi_iter
                        =work.end(); Bisi_iter--;

                    // Keep iterating until Bisi equals the first element in the
                    // work list. This means we have computed B1s1, B2s2, ...,
                    // Bksk.
                    Bisj_iter=work.begin();
                    typename std::list <X_Vector>::const_iterator si_iter
                        =oldS.end(); si_iter--;
                    typename std::list <X_Vector>::const_iterator yi_iter
                        =oldY.end(); yi_iter--;
                    typename std::list <X_Vector>::const_iterator sj_iter
                        =oldS.begin();
                    while(true){

                        // Create some reference to our iterators that are
                        // easier to work with
                        X_Vector const & si=*si_iter;
                        X_Vector const & yi=*yi_iter;
                        X_Vector const & Bisi=*Bisi_iter;

                        // Determine <Bi si,si>
                        Real inner_Bisi_si=X::innr(Bisi,si);

                        // Determine <yi,si>
                        Real inner_yi_si=X::innr(yi,si);

                        // Determine <si,Bi dx>
                        Real inner_si_Bidx=X::innr(si,result);

                        // Determine <yi,dx>
                        Real inner_yi_dx=X::innr(yi,dx);

                        // Determine -<si,Bi dx>/<Bi si,si> Bisi + Bi dx.
                        // Store in Bi dx.  This will become B_{i+1} dx.
                        X::axpy(-inner_si_Bidx/inner_Bisi_si,Bisi,result);

                        // Determine <yi,dx>/<yi,si> yi + w where we calculated
                        // w in the line above.  This completes the calculation
                        // of B_{i+1} dx
                        X::axpy(inner_yi_dx/inner_yi_si,yi,result);

                        // Check whether or not we've calculated B_{i+1} dx for
                        // the last time
                        if(Bisi_iter==work.begin()) break;

                        // Begin the calculation of B_{i+1}sj
                        while(si_iter!=sj_iter){
                            // Add some additional references to the iterators 
                            X_Vector const & sj=*sj_iter;
                            X_Vector & Bisj=*Bisj_iter;

                            // Determine <si,Bisj>
                            Real inner_si_Bisj=X::innr(si,Bisj);

                            // Determine <yi,sj>
                            Real inner_yi_sj=X::innr(yi,sj);

                            // Determine -<si,Bisj>/<Bisi,si> Bisi + Bisj
                            // Store in Bisj.  This will become B_{i+1}sj.
                            X::axpy(-inner_si_Bisj/inner_Bisi_si,Bisi,Bisj);

                            // Determine <yi,sj>/<yi,si> yi + w where we 
                            // calculated w in the line above.  This completes 
                            // the computation of B_{i+1}sj.
                            X::axpy(inner_yi_sj/inner_yi_si,yi,Bisj);

                            // Change j to be j-1 and adjust Bisj and sj 
                            // accordingly
                            sj_iter++;
                            Bisj_iter++;
                        }

                        // At this point, we've computed all Bisj entries on the
                        // current row.  As a result, we increment i and set j 
                        // to be k.  This requires us to modify si, yi, sj, 
                        // Bisj, and Bisi accordingly.
                        
                        // Increment i and adjust si
                        si_iter--;

                        // Increment i and adjust yi
                        yi_iter--;

                        // Set j=k and adjust sj
                        sj_iter=oldS.begin();

                        // Set j=k, increment i, and adjust Bisj
                        Bisj_iter=work.begin();

                        // Increment i and adjust Bisi
                        Bisi_iter--;
                    }
                }
            };

            // The SR1 Hessian approximation.  
            class SR1 : public Operator <Real,XX,XX> {
            private:
                // Messaging device in case the quasi-Newton information is bad
                Messaging const & msg;

                // Stored quasi-Newton information
                std::list<X_Vector> const & oldY;
                std::list<X_Vector> const & oldS;
            public:
                SR1(
                    Messaging const & msg_,
                    typename State::t const & state
                ) : msg(msg_), oldY(state.oldY), oldS(state.oldS) {};
                
                // Operator interface
                void operator () (X_Vector const & dx,X_Vector & result) const {

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the SR1 Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work;
                    for(int i=0;i<oldY.size();i++)
                        work.emplace_back(std::move(X::init(dx)));

                    // If we have no vectors in our history, we return the 
                    // direction
                    X::copy(dx,result);
                    if(oldY.size() == 0) return;

                    // Othwerwise, we copy all of the trial step differences 
                    // into the work space
                    typename std::list <X_Vector>::iterator Bisj_iter
                        =work.begin();
                    typename std::list <X_Vector>::const_iterator sk_iter
                        =oldS.begin();
                    while(Bisj_iter!=work.end())
                        X::copy((*sk_iter++),(*Bisj_iter++));

                    // Keep track of the element Bisi
                    typename std::list <X_Vector>::const_iterator Bisi_iter
                        =work.end(); Bisi_iter--;

                    // Keep iterating until Bisi equals the first element in the
                    // work list. This means we have computed B1s1, B2s2, ...,
                    // Bksk.
                    Bisj_iter=work.begin();
                    typename std::list <X_Vector>::const_iterator si_iter
                        =oldS.end(); si_iter--;
                    typename std::list <X_Vector>::const_iterator yi_iter
                        =oldY.end(); yi_iter--;
                    typename std::list <X_Vector>::const_iterator sj_iter
                        =oldS.begin();
                    while(true){

                        // Create some reference to our iterators that are 
                        // easier to work with
                        X_Vector const & si=*si_iter;
                        X_Vector const & yi=*yi_iter;
                        X_Vector const & Bisi=*Bisi_iter;

                        // Determine <yi,dx>
                        Real inner_yi_dx=X::innr(yi,dx);

                        // Determine <Bisi,dx>
                        Real inner_Bisi_dx=X::innr(Bisi,dx);

                        // Determine <yi,si>
                        Real inner_yi_si=X::innr(yi,si);

                        // Determine <Bisi,si>
                        Real inner_Bisi_si=X::innr(Bisi,si);

                        // Determine (<yi,p>-<Bisi,dx>) / (<y_i,s_i>-<Bisi,si>).
                        // Store in alpha
                        Real alpha=(inner_yi_dx-inner_Bisi_dx)
                            / (inner_yi_si-inner_Bisi_si);

                        // Determine alpha y_i + Bi dx.  Store in result (which
                        // accumulate Bi dx).
                        X::axpy(alpha,yi,result);

                        // Then, add -alpha*Bisi to this result
                        X::axpy(-alpha,Bisi,result);

                        // Check whether or not we've calculated B_{i+1}p for 
                        // the last time
                        if(Bisi_iter==work.begin()) break;

                        // Begin the calculation of B_{i+1}sj
                        while(si_iter!=sj_iter){
                            // Add some additional references to the iterators 
                            X_Vector const & sj=*sj_iter;
                            X_Vector & Bisj=*Bisj_iter;

                            // Determine <yi,sj>
                            Real inner_yi_sj=X::innr(yi,sj);

                            // Determine <Bisi,sj>
                            Real inner_Bisi_sj=X::innr(Bisi,sj);

                            // Determine
                            // (<yi,dx>-<Bisi,dx>) / (<yi,si>-<Bisi,si>).
                            // Store in beta.
                            //
                            // CHECK THIS FORMULA
                            Real beta = (inner_yi_sj-inner_Bisi_sj) /
                                (inner_yi_si-inner_Bisi_si);
                        
                            // Determine beta y_i + Bisj.  Store in Bisj. 
                            X::axpy(beta,yi,Bisj);

                            // Add -beta*Bisi to this result
                            X::axpy(-beta,Bisi,Bisj);

                            // Change j to be j-1 and adjust Bisj and sj 
                            // accordingly
                            sj_iter++;
                            Bisj_iter++;
                        }

                        // At this point, we've computed all Bisj entries on the
                        // current row.  As a result, we increment i and set j 
                        // to be k.  This requires us to modify si, yi, sj, 
                        // Bisj, and Bisi accordingly.
                        
                        // Increment i and adjust si
                        si_iter--;

                        // Increment i and adjust yi
                        yi_iter--;

                        // Set j=k and adjust sj
                        sj_iter=oldS.begin();

                        // Set j=k, increment i, and adjust Bisj
                        Bisj_iter=work.begin();

                        // Increment i and adjust Bisi
                        Bisi_iter--;
                    }
                }
            };

            // The inverse BFGS operator 
            class InvBFGS : public Operator <Real,XX,XX> {
            private:
                // Messaging device in case the quasi-Newton information is bad
                Messaging const & msg;

                // Stored quasi-Newton information
                std::list <X_Vector> const & oldY;
                std::list <X_Vector> const & oldS;
            public:
                InvBFGS(
                    Messaging const & msg_,
                    typename State::t const & state
                ) : msg(msg_), oldY(state.oldY), oldS(state.oldS) {};
                
                // Operator interface
                void operator () (X_Vector const & dx,X_Vector & result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the inverse BFGS operator, the number "
                            "of stored gradient differences must equal the "
                            "number of stored trial step differences.");
                    
                    // As a safety check, insure that the inner product between
                    // all the (s,y) pairs is positive
                    typename std::list <X_Vector>::const_iterator y0
                        = oldY.begin();
                    typename std::list <X_Vector>::const_iterator s0
                        = oldS.begin();
                    while(y0!=oldY.end()){
                        Real inner_y_s=X::innr(*y0++,*s0++);
                        if(inner_y_s <= Real(0.))
                            msg.error("Detected a (s,y) pair in the inverse "
                                "BFGS operator that possesed a nonpositive "
                                "inner product");
                    }

                    // Create two vectors to hold some intermediate calculations
                    std::vector <Real> alpha(oldY.size());
                    std::vector <Real> rho(oldY.size());

                    // Before we begin computing, copy dx to our result 
                    X::copy(dx,result);

                    // In order to compute, we first iterate over all the stored
                    // element in the forward direction.  Then, we iterate over
                    // them backward.
                    typename std::list <X_Vector>::const_iterator y_iter
                        =oldY.begin();
                    typename std::list <X_Vector>::const_iterator s_iter
                        =oldS.begin();
                    Natural i(0);
                    while(y_iter != oldY.end()){
                        // Find y_k, s_k, and their inner product
                        X_Vector const & y_k=*(y_iter++);
                        X_Vector const & s_k=*(s_iter++);
                        rho[i]=Real(1.)/X::innr(y_k,s_k);

                        // Find rho_i <s_i,result>.  Store in alpha_i
                        alpha[i]=rho[i]*X::innr(s_k,result);

                        // result = - alpha_i y_i + result 
                        X::axpy(-alpha[i],y_k,result);

                        // Make sure we don't overwrite alpha and rho
                        i++;
                    }

                    // Assume that H_0 is the identity operator (which may or 
                    // may not work in Hilbert space)

                    // Now, let us iterate backward over our elements to 
                    // complete the computation
                    while(y_iter != oldY.begin()){
                        // Find y_k and s_k
                        X_Vector const & s_k=*(--s_iter);
                        X_Vector const & y_k=*(--y_iter);

                        // beta=rho_i <y_i,result>
                        Real beta= rho[--i] * X::innr(y_k,result);

                        // result=  (alpha_i-beta) s_i + result
                        X::axpy(alpha[i]-beta,s_k,result);
                    }
                }
            };
            
            // The inverse SR1 operator.  In this definition, we take a
            // shortcut and simply use the SR1 Hessian approximation where we
            // swap Y and S.
            class InvSR1 : public Operator <Real,XX,XX> {
            private:
                // Store the SR1 operator
                SR1 sr1;
            public:
                InvSR1(
                    Messaging const & msg,
                    typename State::t const & state
                ) : sr1(msg,state) {};
                void operator () (X_Vector const & dx,X_Vector & result) const{
                    sr1(dx,result);
                }
            };

            // A scalar valued function that overrides the Hessian if need be. 
            struct HessianAdjustedFunction
                : public Optizelle::ScalarValuedFunction <Real,XX>
            {
            private:
                // Hessian approximation
                std::unique_ptr <Operator <Real,XX,XX> > H;

                // Underlying function
                std::unique_ptr <Optizelle::ScalarValuedFunction <Real,XX> > f;

            public:
                // Prevent constructors 
                NO_DEFAULT_COPY_ASSIGNMENT(HessianAdjustedFunction);

                // The constructor determines whether we really need to build
                // a Hessian-vector product or if we use an internal
                // approximation
                HessianAdjustedFunction(
                    Messaging const & msg,
                    typename State::t const & state,
                    typename Functions::t & fns
                ) : H(nullptr), f(std::move(fns.f)) {
                    // Determine the Hessian approximation
                    switch(state.H_type){
                        case Operators::Identity:
                            H.reset(new Identity());
                            break;
                        case Operators::ScaledIdentity:
                            H.reset(new ScaledIdentity (fns,state));
                            break;
                        case Operators::BFGS:
                            H.reset(new BFGS(msg,state));
                            break;
                        case Operators::SR1:
                            H.reset(new SR1(msg,state));
                            break;
                        case Operators::UserDefined:
                            break;
                        default:
                            msg.error("Not a valid Hessian approximation.");
                            break;
                    }
                }

                 // <- f(x) 
                 Real operator () (X_Vector const & x) const {
                    return (*f)(x);
                 }

                 // grad = grad f(x) 
                 void grad(X_Vector const & x,X_Vector & grad) const {
                        f->grad(x,grad);
                 }

                 // H_dx = hess f(x) dx 
                 // This actually computes the Hessian-vector product.  In 
                 // essence, we may want to use a Hessian approximation 
                 // provided by the optimization routines.  The following
                 // routine selects whether or not we use the hessvec 
                 // provided by the user.
                 virtual void hessvec(
                     X_Vector const & x,
                     X_Vector const & dx,
                     X_Vector & H_dx 
                 ) const {
                     if(H.get()!=nullptr) 
                        (*H)(dx,H_dx);
                     else
                        f->hessvec(x,dx,H_dx);
                 }
            };

            // Check that all the functions are defined
            static void check(Messaging const & msg,t const & fns) {
                // Check that objective function exists 
                if(fns.f.get()==nullptr)
                    msg.error("Missing an objective function definition.");
                
                // Check that objective function modifications exists 
                if(fns.f_mod.get()==nullptr)
                    msg.error("Missing an objective function modification "
                        "definition.");
                
                // Check that a preconditioner exists 
                if(fns.PH.get()==nullptr)
                    msg.error("Missing a preconditioner definition.");
            }

            // Initialize any missing functions for just unconstrained
            // optimization.
            static void init_(
                Messaging const & msg,
                typename State::t const & state,
                t & fns
            ) {
                // Create the objective modifications
                fns.f_mod.reset(
                    new ScalarValuedFunctionModifications <Real,XX>());

                // Determine the preconditioner
                switch(state.PH_type){
                    case Operators::Identity:
                        fns.PH.reset(new Identity());
                        break;
                    case Operators::InvBFGS:
                        fns.PH.reset(new InvBFGS(msg,state));
                        break;
                    case Operators::InvSR1:
                        fns.PH.reset(new InvSR1(msg,state));
                        break;
                    case Operators::UserDefined:
                        if(fns.PH.get()==nullptr)
                            msg.error("An externally defined preconditioner "
                                "must be provided explicitly.");
                        break;
                    default:
                        msg.error("Not a valid Hessian approximation.");
                        break;
                }

                // If a trust-region operator has not been provided, use the
                // identity.
                if(fns.TRS.get()==nullptr)
                    fns.TRS.reset(new Identity());

                // Check that all functions are defined (namely, the 
                // objective).
                check(msg,fns);

                // Modify the objective function if necessary
                fns.f.reset(new HessianAdjustedFunction(msg,state,fns));
            }

            // Initialize any missing functions 
            static void init(
                Messaging const & msg,
                typename State::t const & state,
                t & fns
            ) {
                Unconstrained <Real,XX>::Functions::init_(msg,state,fns);
            }
        };

        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
            // Disallow constructors
            NO_CONSTRUCTORS(Printer);

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                LineSearchDirection::t const & dir=state.dir;

                // Basic information
                out.emplace_back(Utility::atos("Iter"));
                out.emplace_back(Utility::atos("f(x)"));
                out.emplace_back(Utility::atos("merit(x)"));
                out.emplace_back(Utility::atos("||grad||"));
                out.emplace_back(Utility::atos("||dx||"));

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    out.emplace_back(Utility::atos("KryIter"));
                    out.emplace_back(Utility::atos("KryErr"));
                    out.emplace_back(Utility::atos("KryWhy"));
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    out.emplace_back(Utility::atos("LSIter"));
                }

                // In case we're using a trust-region method 
                if(algorithm_class==AlgorithmClass::TrustRegion) {
                    out.emplace_back(Utility::atos("ared"));
                    out.emplace_back(Utility::atos("pred"));
                    out.emplace_back(Utility::atos("ared/pred"));
                }
            }

            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
            }

            // Gets the state information for output
            static void getState_(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                bool const & noiter,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                X_Vector const & grad=state.grad;
                Natural const & iter=state.iter;
                Real const & f_x=state.f_x;
                Natural const & krylov_iter=state.krylov_iter;
                Real const & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Natural const & linesearch_iter=state.linesearch_iter;
                Real const & ared=state.ared;
                Real const & pred=state.pred;
                Real const & alpha=state.alpha;
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                LineSearchDirection::t const & dir=state.dir;
                Natural const & rejected_trustregion=state.rejected_trustregion;

                // Figure out if we're at the absolute beginning of the
                // optimization.  We have to be a little saavy about this
                // since we could be on the first iteration, but in the
                // middle of a line-search or trust-region method and
                // still want to output things
                bool opt_begin = (iter==1) &&
                    ((algorithm_class == AlgorithmClass::LineSearch && 
                        linesearch_iter==0) ||
                    ((algorithm_class == AlgorithmClass::TrustRegion ||
                      algorithm_class == AlgorithmClass::UserDefined) && 
                        rejected_trustregion == 0));

                // Determine some extra diagnostic information
                Real merit_x=f_mod.merit(x,f_x);
                Real norm_dx=sqrt(X::innr(dx,dx));
                X_Vector grad_diag(X::init(grad));
                    f_mod.grad_diag(x,grad,grad_diag);
                Real norm_grad=sqrt(X::innr(grad_diag,grad_diag));


                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                if(!noiter)
                    out.emplace_back(Utility::atos(iter));
                else
                    out.emplace_back(Utility::atos("*"));
                out.emplace_back(Utility::atos(f_x));
                out.emplace_back(Utility::atos(merit_x));
                out.emplace_back(Utility::atos(norm_grad));
                if(!opt_begin) {
                    if(algorithm_class==AlgorithmClass::LineSearch)
                        out.emplace_back(Utility::atos(alpha*norm_dx));
                    else
                        out.emplace_back(Utility::atos(norm_dx));
                } else
                    out.emplace_back(Utility::blankSeparator);

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    if(!opt_begin) {
                        out.emplace_back(Utility::atos(krylov_iter));
                        out.emplace_back(Utility::atos(krylov_rel_err));
                        out.emplace_back(Utility::atos(krylov_stop));
                    } else 
                        for(Natural i=0;i<3;i++)
                            out.emplace_back(Utility::blankSeparator);
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    if(!opt_begin)
                        out.emplace_back(Utility::atos(linesearch_iter));
                    else 
                        out.emplace_back(Utility::blankSeparator);
                }
                
                // In case we're using a trust-region method
                if(algorithm_class==AlgorithmClass::TrustRegion) {
                    if(!opt_begin) {
                        out.emplace_back(Utility::atos(ared));
                        out.emplace_back(Utility::atos(pred));
                        out.emplace_back(Utility::atos(ared/pred));
                    } else  
                        for(Natural i=0;i<3;i++)
                            out.emplace_back(Utility::blankSeparator);
                }

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=Utility::blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                bool const & noiter,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
            }

            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                // Create some shortcuts
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                LineSearchDirection::t const & dir=state.dir;

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    out.emplace_back(Utility::atos("KrySubItr"));
                    out.emplace_back(Utility::atos("KryTotItr"));
                    out.emplace_back(Utility::atos("KrySubErr"));
                }
            }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                // Create some shortcuts
                Natural const & krylov_iter=state.krylov_iter;
                Natural const & krylov_iter_total=state.krylov_iter_total;
                Real const & krylov_rel_err=state.krylov_rel_err;

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                out.emplace_back(Utility::atos(krylov_iter));
                out.emplace_back(Utility::atos(krylov_iter_total));
                out.emplace_back(Utility::atos(krylov_rel_err));
                
                // If we needed to do blank insertions, overwrite the elements
                // with nothing
                if(blank) {
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        x->clear();
                }
            }

            // Combines all of the Krylov information
            static void getKrylov(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms);

            // Checks a set of stopping conditions
            static StoppingCondition::t checkStop(
                typename Functions::t const & fns, 
                typename State::t const & state
            ){
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & dx=state.dx;
                Real const & norm_gradtyp=state.norm_gradtyp;
                Real const & norm_dxtyp=state.norm_dxtyp;
                Natural const & iter=state.iter;
                Natural const & iter_max=state.iter_max;
                Real const & eps_grad=state.eps_grad;
                Real const & eps_dx=state.eps_dx;

                // Find both the norm of the gradient and the step
                X_Vector grad_stop(X::init(grad));
                f_mod.grad_stop(x,grad,grad_stop);
                const Real norm_grad=sqrt(X::innr(grad_stop,grad_stop));
                const Real norm_dx=sqrt(X::innr(dx,dx));

                // Check if we've exceeded the number of iterations
                if(iter>=iter_max)
                    return StoppingCondition::MaxItersExceeded;

                // Check whether the change in the step length has become too
                // small relative to some typical step
                if(norm_dx < eps_dx*norm_dxtyp)
                    return StoppingCondition::RelativeStepSmall;
                
                // Check whether the norm is small relative to some typical
                // gradient
                if(norm_grad < eps_grad*norm_gradtyp)
                    return StoppingCondition::RelativeGradientSmall;

                // Otherwise, return that we're not converged 
                return StoppingCondition::NotConverged;
            }

            // Sets up the Hessian operator for use in the Krylov methods.  In
            // other words, this sets up the application H(x)dx.
            struct HessianOperator : public Operator <Real,XX,XX> {
            private:
                // Store the objective
                ScalarValuedFunction <Real,XX> const & f;

                // Objective modifications
                ScalarValuedFunctionModifications <Real,XX> const & f_mod;

                // Store a reference to the base of the Hessian-vector product
                X_Vector const & x;

                // Allocate memory for the Hessian modification
                mutable X_Vector H_dx;

            public:
                // Take in the objective and the base point during construction 
                HessianOperator(
                    ScalarValuedFunction <Real,XX> const & f_,
                    ScalarValuedFunctionModifications <Real,XX> const & f_mod_,
                    X_Vector const & x_)
                : f(f_), f_mod(f_mod_), x(x_), H_dx(X::init(x_))
                {}

                // Basic application
                void operator () (X_Vector const & dx,X_Vector & Hdx_step)
                    const
                {
                    f.hessvec(x,dx,H_dx);
                    f_mod.hessvec_step(x,dx,H_dx,Hdx_step);
                }
            };
        
            // Checks whether we accept or reject a step
            static bool checkStep(
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                X_Vector const & grad=state.grad;
                Real const & eta1=state.eta1;
                Real const & eta2=state.eta2;
                Real const & f_x=state.f_x;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Real & delta=state.delta;
                Real & ared=state.ared;
                Real & pred=state.pred;
                Real & f_xpdx=state.f_xpdx;
                
                // Allocate memory for temporaries that we need
                X_Vector x_p_dx(X::init(x));

                // Determine x+dx 
                X::copy(dx,x_p_dx);
                X::axpy(Real(1.),x,x_p_dx);

                // Determine merit(x)
                Real merit_x = f_mod.merit(x,f_x);
                
                // Determine H(x)dx
                X_Vector H_dx(X::init(x));
                    f.hessvec(x,dx,H_dx);
                X_Vector Hdx_step(X::init(x));
                    f_mod.hessvec_step(x,dx,H_dx,Hdx_step);

                // Determine the gradient
                X_Vector grad_step(X::init(x));
                    f_mod.grad_step(x,grad,grad_step);

                // Calculate the model,
                // m(dx) = f(x) + < grad,dx > + < H(x)dx,dx >
                Real model_dx= merit_x + X::innr(grad_step,dx)
                    + Real(.5)*X::innr(Hdx_step,dx);

                // Determine the merit function evaluated at x+dx
                f_xpdx=f(x_p_dx);
                Real merit_xpdx=f_mod.merit(x_p_dx,f_xpdx);

                // Determine the norm of the step
                Real norm_dx = sqrt(X::innr(dx,dx));

                // Determine the reductions
                ared = merit_x - merit_xpdx;
                pred = merit_x - model_dx;

                // Add a safety check in case we don't actually minimize the TR
                // subproblem correctly. This could happen for a variety of
                // reasons.  Most notably, if we do not correctly calculate the
                // Hessian approximation, we could have a nonsymmetric 
                // approximation.  In that case, truncated-CG will exit, but 
                // has an undefined result.  In the case that the actual 
                // reduction also increases, rho could have an extraneous 
                // positive value.  Hence, we require an extra check.
                if(model_dx > merit_x){
                    delta = norm_dx/Real(2.);
                    return false;
                }

                // Update the trust region radius and return whether or not we
                // accept the step
                if(ared >= eta2*pred){
                    // Increase the size of the trust-region if the Krylov
                    // solver reached the boundary.
                    if( krylov_stop==KrylovStop::NegativeCurvature ||
                        krylov_stop==KrylovStop::TrustRegionViolated
                    ) 
                        delta*=Real(2.);
                    return true;
                } else if(ared >= eta1*pred && ared < eta2*pred)
                    return true;
                else {
                    delta = norm_dx/Real(2.);
                    return false;
                }
            }
        
            // Finds the trust-region step
            static void getStepTR(
                Messaging const & msg,
                StateManipulator <Unconstrained <Real,XX> > const & smanip,
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                Operator <Real,XX,XX> const & TRS=*(fns.TRS);
                Real const & eps_dx=state.eps_dx;
                Real const & eps_krylov=state.eps_krylov;
                Natural const & krylov_iter_max=state.krylov_iter_max;
                Natural const & krylov_orthog_max=state.krylov_orthog_max;
                Real const & delta=state.delta;
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                Real const & norm_dxtyp=state.norm_dxtyp;
                KrylovSolverTruncated::t const & krylov_solver
                    = state.krylov_solver;
                Natural & rejected_trustregion=state.rejected_trustregion;
                X_Vector & dx=state.dx;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                Real & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                Natural & history_reset=state.history_reset;
                
                // Allocate some memory for the scaled trial step and the
                // trust-region center
                X_Vector x_tmp1(X::init(x));
                X_Vector dx_cp(X::init(x));
                X_Vector grad_step(X::init(x));
                X_Vector minus_grad(X::init(x));

                // Find -grad f(x) 
                f_mod.grad_step(x,grad,grad_step);
                X::copy(grad_step,minus_grad);
                X::scal(Real(-1.),minus_grad);

                // Continue to look for a step until one comes back as valid
                for(rejected_trustregion=0;
                    true; 
                ) {
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeGetStep);

                    // Use truncated the truncated Krylov solver to find a 
                    // new trial step
                    HessianOperator H(f,f_mod,x);

                    // Set the trust-region center
                    X::zero(x_tmp1);

                    // Keep track of the residual errors
                    Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                    Real residual_err(std::numeric_limits <Real>::quiet_NaN());

                    switch(krylov_solver) {
                    // Truncated conjugate direction
                    case KrylovSolverTruncated::ConjugateDirection:
                        truncated_cd(
                            H,
                            minus_grad,
                            PH,
                            TRS,
                            eps_krylov,
                            krylov_iter_max,
                            krylov_orthog_max,
                            delta,
                            x_tmp1,
			    false,
                            dx,
                            dx_cp,
                            residual_err0,
                            residual_err,
                            krylov_iter,
                            krylov_stop);
                        break;

                    // Truncated MINRES 
                    case KrylovSolverTruncated::MINRES:
                        truncated_minres(
                            H,
                            minus_grad,
                            PH,
                            TRS,
                            eps_krylov,
                            krylov_iter_max,
                            krylov_orthog_max,
                            delta,
                            x_tmp1,
                            dx,
                            dx_cp,
                            residual_err0,
                            residual_err,
                            krylov_iter,
                            krylov_stop);

                        // Force a descent direction
                        if(X::innr(dx,grad) > 0) X::scal(Real(-1.),dx);
                        if(X::innr(dx_cp,grad) > 0) X::scal(Real(-1.),dx_cp);
                        break;
                    }
                    krylov_rel_err = residual_err 
                        / (std::numeric_limits <Real>::epsilon()+residual_err0);
                    krylov_iter_total += krylov_iter;

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::BeforeActualVersusPredicted);

                    // Check whether the step is good
                    if(checkStep(fns,state))
                        break;
                    else
                        rejected_trustregion++;
                    
                    // If the number of rejected steps is above the
                    // history_reset threshold, destroy the quasi-Newton
                    // information
                    if(rejected_trustregion > history_reset){
                        oldY.clear();
                        oldS.clear();
                    }

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::AfterRejectedTrustRegion);

                    // Alternatively, check if the step becomes so small
                    // that we're not making progress.  In this case, break
                    // and allow the stopping conditions to terminate
                    // optimization.  We use a zero length step so that we
                    // do not modify the current iterate.
                    Real norm_dx = sqrt(X::innr(dx,dx));
                    if(norm_dx < eps_dx*norm_dxtyp) {
                        X::scal(Real(0.),dx);
                        norm_dx=Real(0.);
                        break;
                    }
                } 
            }
        
            // Steepest descent search direction
            static void SteepestDescent(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector & dx=state.dx;

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);

                // We take the steepest descent direction and apply the
                // preconditioner.
                PH(grad_step,dx);
                X::scal(Real(-1.),dx);
            }
    
            // Nonlinear Conjugate Gradient
            static void NonlinearCG(
                typename NonlinearCGDirections::t const & dir,
                typename Functions::t const & fns,
                typename State::t & state
            ) {
            
                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                Real const & alpha=state.alpha;
                X_Vector & dx_old=state.dx_old;
                Natural & iter=state.iter;
                X_Vector & dx=state.dx;

                // Scale dx by 1/alpha.  In our algorithms, we always stored
                // alpha*dx in order to better integrate with trust-region
                // algorithms.  In nonlinear-CG, the formulas assume that
                // we have not stored the scaled direction.  Hence, we scale
                // it here and fix it at the end of the routine.
                X::scal(1./alpha,dx_old);

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);

                // If we're on the first iterations, we take the steepest
                // descent direction
                if(iter==1) SteepestDescent(fns,state);

                // On subsequent iterations, we take the specified direction
                else {
                    // Find the momentum parameter
                    Real beta(std::numeric_limits<Real>::quiet_NaN());
                    switch(dir) {
                    case NonlinearCGDirections::FletcherReeves:
                        beta=FletcherReeves(fns,state);
                        break;
                    case NonlinearCGDirections::PolakRibiere:
                        beta=PolakRibiere(fns,state);
                        break;
                    case NonlinearCGDirections::HestenesStiefel:
                        beta=HestenesStiefel(fns,state);
                        break;
                    }

                    // Find -PH grad+beta*dx_old.  
                    PH(grad_step,dx);
                    X::scal(Real(-1.),dx);
                    X::axpy(beta,dx_old,dx);

                    // We don't ever check the strong-Wolfe conditions, so
                    // hard check that we have a descent direction
                    if(X::innr(dx,grad_step) > 0) X::scal(Real(-1.),dx);
                }

                // Undo the scaling of the previous search direction
                X::scal(alpha,dx_old);
            }

            // Fletcher-Reeves CG search direction
            static Real FletcherReeves(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & grad_old=state.grad_old;

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector grad_old_step(X::init(grad));
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Apply the preconditioner to the gradients 
                X_Vector PH_grad_step(X::init(grad_step));
                    PH(grad_step,PH_grad_step);
                X_Vector PH_grad_old_step(X::init(grad_old_step));
                    PH(grad_old_step,PH_grad_old_step);

                // Return the momentum parameter
                return X::innr(grad_step,PH_grad_step)
                    / X::innr(grad_old_step,PH_grad_old_step);
            }
        
            // Polak-Ribiere CG search direction
            static Real PolakRibiere(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & grad_old=state.grad_old;

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector grad_old_step(X::init(grad));
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Find grad-grad_old 
                X_Vector grad_m_gradold(X::init(grad));
                X::copy(grad_step,grad_m_gradold);
                X::axpy(Real(-1.),grad_old_step,grad_m_gradold);
                
                // Apply the preconditioner to the gradients 
                X_Vector PH_grad_step(X::init(grad_step));
                    PH(grad_step,PH_grad_step);
                X_Vector PH_grad_old_step(X::init(grad_old_step));
                    PH(grad_old_step,PH_grad_old_step);
                    
                // Return the momentum parameter
                return X::innr(PH_grad_step,grad_m_gradold)
                    / X::innr(PH_grad_old_step,grad_old_step);
            }
            
            // Hestenes-Stiefel search direction
            static Real HestenesStiefel(
                typename Functions::t const & fns,
                typename State::t & state
            ) {

                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & grad_old=state.grad_old;
                X_Vector const & dx_old=state.dx_old;

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector grad_old_step(X::init(grad));
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Find grad-grad_old 
                X_Vector grad_m_gradold(X::init(grad));
                X::copy(grad_step,grad_m_gradold);
                X::axpy(Real(-1.),grad_old_step,grad_m_gradold);
                
                // Apply the preconditioner to the gradient
                X_Vector PH_grad_step(X::init(grad_step));
                    PH(grad_step,PH_grad_step);
                    
                // Return the momentum parameter.
                Real beta=X::innr(PH_grad_step,grad_m_gradold)
                    / X::innr(dx_old,grad_m_gradold);
                return beta < Real(0.) ? Real(0.) : beta;
            }

            // BFGS search direction
            static void BFGS(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                
                // Create some shortcuts 
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector & dx=state.dx;

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);

                // Create the inverse BFGS operator
                typename Functions::InvBFGS Hinv(msg,state); 

                // Apply the inverse BFGS operator to the gradient
                Hinv(grad_step,dx);

                // Negate the result
                X::scal(Real(-1.),dx);
            }

            // Compute a Golden-Section search between 0 and alpha0. 
            static typename LineSearchTermination::t goldenSection(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                Natural const & iter_max=state.linesearch_iter_max;
                Real const & alpha0=state.alpha0;
                Natural & iter_total=state.linesearch_iter_total;
                Natural & iter=state.linesearch_iter;
                Real & f_xpdx=state.f_xpdx;
                Real & alpha=state.alpha;
                
                // Create one work element that holds x+mu dx or x+lambda dx 
                X_Vector x_p_dx(X::init(x));

                // Find 1 over the golden ratio
                Real beta=Real(2./(1.+sqrt(5.)));

                // Find a bracket for the linesearch such that a < b
                Real a=Real(0.);
                Real b=alpha0;

                // Find two new points between a and b, mu and lambda,
                // such that lambda < mu
                Real lambda=a+(1.-beta)*(b-a);
                Real mu=a+beta*(b-a);

                // Find the merit value at mu and labmda 

                // mu 
                X::copy(x,x_p_dx);
                X::axpy(mu,dx,x_p_dx);
                Real f_mu=f(x_p_dx);
                Real merit_mu=f_mod.merit(x_p_dx,f_mu);

                // lambda
                X::copy(x,x_p_dx);
                X::axpy(lambda,dx,x_p_dx);
                Real f_lambda=f(x_p_dx);
                Real merit_lambda=f_mod.merit(x_p_dx,f_lambda);

                // Search for a fixed number of iterations.  Note, since we
                // already evaluated the objective twice above, at mu and
                // lambda, we've already done two iterations.  Hence, we assume
                // that iter_max >=2.
                for(iter=2;iter<iter_max;iter++){

                    // If the merit is greater on the left, bracket on the
                    // right.  Alternatively, it's possible that we're going to
                    // generate a NaN on the right.  This means that
                    // merit_mu=NaN.  In this case we want to bracket on the
                    // left.  Since merit_lambda > merit_mu will return false 
                    // when merit_mu is a NaN, we should be safe.
                    if(merit_lambda > merit_mu){
                        a=lambda;
                        lambda=mu;
                        merit_lambda=merit_mu;
                        mu=a+beta*(b-a);

                        X::copy(x,x_p_dx);
                        X::axpy(mu,dx,x_p_dx);
                        f_mu=f(x_p_dx);
                        merit_mu=f_mod.merit(x_p_dx,f_mu);

                    // Otherwise, the objective is greater on the right, so
                    // bracket on the left
                    } else {
                        b=mu;
                        mu=lambda;
                        merit_mu=merit_lambda;
                        lambda=a+(1-beta)*(b-a);
                
                        X::copy(x,x_p_dx);
                        X::axpy(lambda,dx,x_p_dx);
                        f_lambda=f(x_p_dx);
                        merit_lambda=f_mod.merit(x_p_dx,f_lambda);
                    }
                }

                // Keep track of the total number of iterations
                iter_total += iter;

                // Once we're finished narrowing in on a solution, take our best
                // guess for the line search parameter and adjust the step
                // length
                alpha=merit_lambda < merit_mu ? lambda : mu;

                // Save the objective value at this step
                f_xpdx=merit_lambda < merit_mu ? f_lambda : f_mu;

                // Determine whether or not we hit the bounds
                if(a==Real(0.))
                    return LineSearchTermination::Min;
                else if(b==alpha0)
                    return LineSearchTermination::Max;
                else
                    return LineSearchTermination::Between;
            }
            
            // This doesn't really do anything save setting the line-search
            // parameter to the be the base line-search parameter and evaluating
            // the objective at x+alpha dx.  Really, we're using the safe
            // guard procedure that checks the sufficient decrease condition
            // in order to do the line-search.
            static void backTracking(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                Real const & alpha0=state.alpha0;
                Natural & iter_total=state.linesearch_iter_total;
                Natural & iter=state.linesearch_iter;
                Real & f_xpdx=state.f_xpdx;
                Real & alpha=state.alpha;
                            
                // Set alpha to the base alpha 
                alpha=alpha0;
               
                // Determine x+alpha dx 
                X_Vector x_p_adx(X::init(x));
                    X::copy(x,x_p_adx);
                    X::axpy(alpha,dx,x_p_adx);
    
                // Determine the objective function evaluated at x+dx
                f_xpdx=f(x_p_adx);

                // Set the number of line-search iterations
                iter=1;
                iter_total+=iter;
            }

            // Find the line search parameter based on the 2-point approximation
            // from Barzilai and Borwein
            static void twoPoint(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & x_old=state.x_old;
                X_Vector const & grad_old=state.grad_old;
                LineSearchKind::t const & kind=state.kind;
                Real & alpha=state.alpha;
                X_Vector & dx=state.dx;
                Natural & iter_total=state.linesearch_iter_total;
                Natural & iter=state.linesearch_iter;
                Real & f_xpdx=state.f_xpdx;

                // Find delta_x
                X_Vector delta_x(X::init(x));
                    X::copy(x,delta_x);
                    X::axpy(Real(-1.),x_old,delta_x);

                // Determine the gradient for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
                
                X_Vector grad_old_step(X::init(grad));
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Find delta_grad
                X_Vector delta_grad(X::init(x));
                    X::copy(grad_step,delta_grad);
                    X::axpy(Real(-1.),grad_old_step,delta_grad);

                // Find alpha
                if(kind==LineSearchKind::TwoPointA)
                    alpha=X::innr(delta_x,delta_grad) 
                        / X::innr(delta_grad,delta_grad);
                else if(kind==LineSearchKind::TwoPointB)
                    alpha=X::innr(delta_x,delta_x)/X::innr(delta_x,delta_grad);

                // Save the objective value at this step
                X_Vector x_p_adx(X::init(x));
                    X::copy(x,x_p_adx);
                    X::axpy(alpha,dx,x_p_adx);
                f_xpdx=f(x_p_adx);

                // Since we do one function evaluation, increase the linesearch
                // iteration by one
                iter=1; iter_total++;
            }
            
            // Finds a trial step using a line-search for globalization
            static void getStepLS(
                Messaging const & msg,
                StateManipulator <Unconstrained <Real,XX> > const & smanip,
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                Operator <Real,XX,XX> const & TRS=*(fns.TRS);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                LineSearchDirection::t const & dir=state.dir;
                LineSearchKind::t const & kind=state.kind;
                Natural const & iter=state.iter;
                Real const & f_x=state.f_x;
                Real const & eps_dx=state.eps_dx;
                Real const & norm_dxtyp=state.norm_dxtyp;
                Real const & eps_krylov=state.eps_krylov;
                Natural const & krylov_iter_max=state.krylov_iter_max;
                Natural const & krylov_orthog_max=state.krylov_orthog_max;
                KrylovSolverTruncated::t const & krylov_solver
                    = state.krylov_solver;
                Real const & c1=state.c1;
                X_Vector & dx=state.dx;
                Real & f_xpdx=state.f_xpdx;
                Real & alpha0=state.alpha0;
                Real & alpha=state.alpha;
                Real & krylov_rel_err=state.krylov_rel_err;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                
                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeGetStep);

                // Create the trust-region center 
                X_Vector x_cntr(X::init(x));
                X::zero(x_cntr);

                // Find the line-search direction
                switch(dir){
                case LineSearchDirection::SteepestDescent:
                    SteepestDescent(fns,state);
                    break;
                case LineSearchDirection::FletcherReeves:
                    NonlinearCG(NonlinearCGDirections::FletcherReeves,
                        fns,state);
                    break;
                case LineSearchDirection::PolakRibiere:
                    NonlinearCG(NonlinearCGDirections::PolakRibiere,fns,state);
                    break;
                case LineSearchDirection::HestenesStiefel:
                    NonlinearCG(NonlinearCGDirections::HestenesStiefel,
                        fns,state);
                    break;
                case LineSearchDirection::BFGS:
                    BFGS(msg,fns,state);
                    break;
                case LineSearchDirection::NewtonCG: {
                    HessianOperator H(f,f_mod,x);
                    X_Vector dx_cp(X::init(x));
                    X_Vector grad_step(X::init(grad));
                        f_mod.grad_step(x,grad,grad_step);
                    X_Vector minus_grad(X::init(x));
                        X::copy(grad_step,minus_grad);
                        X::scal(Real(-1.),minus_grad);

                    // Keep track of the residual errors
                    Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                    Real residual_err(std::numeric_limits <Real>::quiet_NaN());

                    switch(krylov_solver) {
                    // Truncated conjugate direction
                    case KrylovSolverTruncated::ConjugateDirection:
                        truncated_cd(
                            H,
                            minus_grad,
                            PH,
                            TRS,
                            eps_krylov,
                            krylov_iter_max,
                            krylov_orthog_max,
                            std::numeric_limits <Real>::infinity(),
                            x_cntr,
                            false,
                            dx,
                            dx_cp,
                            residual_err0,
                            residual_err,
                            krylov_iter,
                            krylov_stop);
                        break;

                    // Truncated MINRES 
                    case KrylovSolverTruncated::MINRES:
                        truncated_minres(
                            H,
                            minus_grad,
                            PH,
                            TRS,
                            eps_krylov,
                            krylov_iter_max,
                            krylov_orthog_max,
                            std::numeric_limits <Real>::infinity(),
                            x_cntr,
                            dx,
                            dx_cp,
                            residual_err0,
                            residual_err,
                            krylov_iter,
                            krylov_stop);

                        // Force a descent direction
                        if(X::innr(dx,grad_step) > 0) X::scal(Real(-1.),dx);
                        if(X::innr(dx_cp,grad_step)>0) X::scal(Real(-1.),dx_cp);
                        break;
                    }
                    krylov_rel_err = residual_err 
                        / (std::numeric_limits <Real>::epsilon()+residual_err0);
                    krylov_iter_total += krylov_iter;
                    break;
                }}

                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeLineSearch);

                // Do the sufficient decrease line-search
                if(LineSearchKind::is_sufficient_decrease(kind) || iter==1) {
                    // Determine merit(x)
                    Real merit_x = f_mod.merit(x,f_x);
                    
                    // Determine the gradient at x
                    X_Vector grad_step(X::init(x));
                        f_mod.grad_step(x,grad,grad_step);
                
                    // Allocate memory for x+alpha dx 
                    X_Vector x_p_adx(X::init(x));

                    // Keep track of whether or not we hit a bound with the
                    // line-search
                    typename LineSearchTermination::t ls_why
                        = LineSearchTermination::Between;

                    // Save the original line search base
                    Real alpha0_orig=alpha0;

                    // Keep track of whether the sufficient decrease condition
                    // is satisfied.
                    bool sufficient_decrease=false;
                    do {
                        // Do the line-search
                        if(kind==LineSearchKind::GoldenSection ||
                            (!LineSearchKind::is_sufficient_decrease(kind) &&
                            iter==1)
                        )
                            ls_why=goldenSection(fns,state);
                        else if(kind==LineSearchKind::BackTracking)
                            backTracking(fns,state);
                        else if(kind==LineSearchKind::Brents) 
                            msg.error("Brent's linesearch is not currently "
                                "implemented.");

                        // Determine x+dx 
                        X::copy(x,x_p_adx);
                        X::axpy(alpha,dx,x_p_adx);
                
                        // Determine the merit function evaluated at x+dx.  This
                        // assumes that the line-search algorithms already
                        // evaluated f(x+alpha dx) and stored it in f_xpdx.
                        Real merit_xpdx=f_mod.merit(x_p_adx,f_xpdx);

                        // Determine if we've satisfied the sufficient decrease
                        // condition.  Also make sure that we don't generate
                        // a NaN
                        sufficient_decrease = 
                            (merit_xpdx==merit_xpdx) &&
                            merit_xpdx < merit_x+c1*alpha*X::innr(grad_step,dx);

                        // If we've not satisfied the sufficient decrease
                        // condition, cut the step
                        if( !sufficient_decrease ) {

                            // Decrease the size of the base line-search
                            // parameter 
                            alpha0/=Real(2.);

                            // Check if the step becomes so small that we're not
                            // making progress.  In this case, take a zero step 
                            // and allow the stopping conditions to exit
                            if(alpha*sqrt(X::innr(dx,dx)) < eps_dx*norm_dxtyp) {
                                X::scal(Real(0.),dx);
                                break;
                            }

                            // Manipulate the state if required
                            smanip(fns,state,
                                OptimizationLocation::AfterRejectedLineSearch);
                        }

                    // Continue as long as we haven't satisfied the sufficient
                    // decrease condition.
                    } while( !sufficient_decrease );
                
                    // If the line-search hit one of the bounds, change the base
                    // line-search parameter.
                    switch(ls_why){
                    case LineSearchTermination::Min:
                        alpha0/=Real(2.);
                        break;
                    case LineSearchTermination::Max:
                        alpha0*=Real(2.);
                        break;
                    case LineSearchTermination::Between:
                        break;
                    }

                    // If we're doing a backtracking line-search, restore the
                    // original line-search parameter base.  Basically, this
                    // line-search doesn't have a mechanism for changing the
                    // base dynamically, so we're assuming the user set a 
                    // reasonable value.  Since the safe guarding modifies
                    // the base line-search parameter, we need to restore it
                    // here.
                    if(kind==LineSearchKind::BackTracking)
                        alpha0=alpha0_orig;

                // Do the line-searches that are not based on sufficient
                // decrease
                } else 
                    twoPoint(fns,state);

                // Adjust the size of the step (apply the line-search 
                // parameter.)
                X::scal(alpha,dx);
            }

            // Finds a new trial step
            static void getStep(
                Messaging const & msg,
                StateManipulator <Unconstrained <Real,XX> > const & smanip,
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;

                // Choose whether we use a line-search or trust-region method
                switch(algorithm_class){
                case AlgorithmClass::TrustRegion:
                    getStepTR(msg,smanip,fns,state);
                    break;
                case AlgorithmClass::LineSearch:
                    getStepLS(msg,smanip,fns,state);
                    break;
                case AlgorithmClass::UserDefined:
                    smanip(fns,state,OptimizationLocation::GetStep);
                    break;
                }
            }

            // Updates the quasi-Newton information
            static void updateQuasi(
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Exit immediately if we're not using a quasi-Newton method
                if(state.stored_history==0) return;

                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & x_old=state.x_old;
                X_Vector const & grad_old=state.grad_old;
                const Operators::t& PH_type=state.PH_type;
                const Operators::t& H_type=state.H_type;
                LineSearchDirection::t const & dir=state.dir;
                std::list <X_Vector>& oldY=state.oldY;
                std::list <X_Vector>& oldS=state.oldS;
               
                // Allocate some temp storage for y and s
                X_Vector s(X::init(x));
                X_Vector y(X::init(x));

                // Find s = x-x_old
                X::copy(x,s);
                X::axpy(Real(-1.),x_old,s);
                
                // Determine the gradient for the quasi-Newton computation 
                X_Vector grad_quasi(X::init(grad));
                    f_mod.grad_quasi(x,grad,grad_quasi);
                X_Vector grad_old_quasi(X::init(grad_old));
                    f_mod.grad_quasi(x,grad_old,grad_old_quasi);

                // Find y = grad - grad_old
                X::copy(grad_quasi,y);
                X::axpy(Real(-1.),grad_old_quasi,y);

                // If we're using BFGS, check that <y,s> > 0
                if((PH_type==Operators::InvBFGS ||
                    H_type==Operators::BFGS ||
                    dir==LineSearchDirection::BFGS)
                    && X::innr(y,s) <= Real(0.))
                    return;

                // Insert these into the quasi-Newton storage
                oldS.emplace_front(std::move(s));
                oldY.emplace_front(std::move(y));

                // Determine if we need to free some memory
                if(oldS.size()>state.stored_history){
                    oldS.pop_back();
                    oldY.pop_back();
                }
            }

            // Solves an optimization problem
            static void getMin_(
                Messaging const & msg,
                StateManipulator <Unconstrained <Real,XX> > const & smanip,
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const &
                    f_mod=*(fns.f_mod);
                X_Vector & x=state.x;
                X_Vector & grad=state.grad;
                X_Vector & dx=state.dx;
                X_Vector & x_old=state.x_old;
                X_Vector & grad_old=state.grad_old;
                X_Vector & dx_old=state.dx_old;
                Real & f_x=state.f_x;
                Real & f_xpdx=state.f_xpdx;
                Real & norm_gradtyp=state.norm_gradtyp;
                Real & norm_dxtyp=state.norm_dxtyp;
                Natural & iter=state.iter;
                StoppingCondition::t & opt_stop=state.opt_stop;
                
                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeginningOfOptimization);

                // Evaluate the objective function and gradient if we've not
                // done so already
                if(f_x != f_x) {

                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation
                        ::BeforeInitialFuncAndGrad);

                    // Sometimes, we can calculate the gradient and objective
                    // simultaneously.  Hence, it's best to calculate the
                    // gradient first and then possibly cache the objective
                    f.grad(x,grad);
                    f_x=f(x);
                    X_Vector grad_stop(X::init(grad));
                        f_mod.grad_stop(x,grad,grad_stop);
                    norm_gradtyp=sqrt(X::innr(grad_stop,grad_stop));

                    // This one is a little funny.  Sometimes, we run into
                    // trouble trying to calculate the initial step.  In this
                    // case, we don't know when to exit due to the relative
                    // step size being small since we've not calculated the
                    // typical step yet.  This safeguards this case.  Basically,
                    // it initially sets the norm of a typical step to be the
                    // norm of the gradient, which is akin to taking a
                    // steepest descent step without globalization.
                    norm_dxtyp=norm_gradtyp;
                
                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::AfterInitialFuncAndGrad);
                }

                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeOptimizationLoop);

                // Primary optimization loop
                do{
                    // Get a new optimization iterate.  
                    getStep(msg,smanip,fns,state);

                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeSaveOld);

                    // Save the old variable, gradient, and trial step.  This
                    // is useful for both CG and quasi-Newton methods.
                    X::copy(x,x_old);
                    X::copy(grad,grad_old);
                    X::copy(dx,dx_old);

                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeStep);

                    // Move to the new iterate
                    X::axpy(Real(1.),dx,x);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::AfterStepBeforeGradient);

                    // Find the new objective value and gradient
                    f_x=f_xpdx;
                    f.grad(x,grad);
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::AfterGradient);
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeQuasi);

                    // Update the quasi-Newton information
                    updateQuasi(fns,state);
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::AfterQuasi);

                    // Increase the iteration
                    iter++;
                    
                    // Check the stopping condition
                    opt_stop=checkStop(fns,state);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::EndOfOptimizationIteration);
                } while(opt_stop==StoppingCondition::NotConverged);
                        
                // Manipulate the state one final time if required
                smanip(fns,state,OptimizationLocation::EndOfOptimization);
            }
            
            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                StateManipulator <Unconstrained <Real,XX> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                StateManipulator <Unconstrained <Real,XX> > const & smanip,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);

                // Check the inputs to the optimization
                State::check(msg,state);

                // Add the output to the state manipulator
                DiagnosticManipulator <Unconstrained<Real,XX> >
                    iomanip(smanip,msg);

                // Minimize the problem
                getMin_(msg,iomanip,fns,state);
            }
        };
    };
    
    // Routines that manipulate and support problems of the form
    // 
    // min_{x \in X} f(x) st g(x) = 0
    //
    // where f : X -> R and g : X -> Y
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class YY
    > 
    struct EqualityConstrained {
        // Disallow constructors
        NO_CONSTRUCTORS(EqualityConstrained);

        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector;

        // This defines a product space between X and Y
        template <typename Real_>
        struct XXxYY {
            typedef std::pair <X_Vector,Y_Vector> Vector;

            // Memory allocation and size setting
            static Vector init(Vector const & x) {
                return std::move(std::pair <X_Vector,Y_Vector> (
                        X::init(x.first),
                        Y::init(x.second)));
            }

            // y <- x (Shallow.  No memory allocation.)
            static void copy(Vector const & x, Vector & y) {
                X::copy(x.first,y.first);
                Y::copy(x.second,y.second);
            }

            // x <- alpha * x
            static void scal(Real_ const & alpha, Vector & x) {
                X::scal(alpha,x.first);
                Y::scal(alpha,x.second);
            }

            // x <- 0 
            static void zero(Vector & x) {
                X::zero(x.first);
                Y::zero(x.second);
            }

            // y <- alpha * x + y
            static void axpy(Real_ const & alpha, Vector const & x, Vector & y){
                X::axpy(alpha,x.first,y.first);
                Y::axpy(alpha,x.second,y.second);
            }

            // innr <- <x,y>
            static Real_ innr(Vector const & x,Vector const & y) {
                return X::innr(x.first,y.first) + Y::innr(x.second,y.second);
            }
        };
        typedef XXxYY <Real> XxY;
        typedef typename XxY::Vector XxY_Vector;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State);

            // Internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t);

                // Equality multiplier (dual variable or Lagrange multiplier)
                Y_Vector y;

                // Step in the equality multiplier 
                Y_Vector dy;

                // The fraction of the total trust-region used for the
                // quasi-norm step
                Real zeta;

                // Trust-region parameter that bounds the error in the
                // predicted reduction
                Real eta0;

                // Penalty parameter for the augmented-Lagrangian
                Real rho;
                
                // Penalty parameter from the last iteration 
                Real rho_old;

                // Fixed increase in the penalty parameter
                Real rho_bar;

                // Stopping tolerance for the norm of the constraints
                Real eps_constr;

                // Inexactness tolerances
                Real xi_qn;     // Quasi-Newton step
                Real xi_pg;     // Projection of the gradient
                Real xi_proj;   // Null-space projection
                Real xi_tang;   // Tangential step
                Real xi_lmh;    // Equality multiplier

                // Sets all the inexactness tolerances 
                void xi_all(Real const & xi) {
                    xi_qn = xi;
                    xi_pg = xi;
                    xi_proj = xi;
                    xi_tang = xi;
                    xi_lmh = xi;
                }

                // Absolute tolerance on the residual of the equality 
                // multiplier solve
                Real xi_lmg;

                // Tolerance for how much error is acceptable after computing
                // the tangential step given the result from the tangential
                // subproblem
                Real xi_4;

                // Residual term in the predicted reduction
                Real rpred;

                // Left preconditioner for the augmented system
                Operators::t PSchur_left_type;
                
                // Right preconditioner for the augmented system
                Operators::t PSchur_right_type;

                // Maximum number of iterations used when solving the augmented
                // system
                Natural augsys_iter_max;

                // How often we restart the augmented system solve
                Natural augsys_rst_freq;
                
                // Equality constraint evaluated at x.  This is used in the
                // quasinormal step as well as in the computation of the
                // linear Taylor series at x in the direciton dx_n.
                Y_Vector g_x;

                // A typical norm for norm_gx.  Generally, we just take
                // the value at the first iteration.
                Real norm_gxtyp;
                
                // Linear Taylor series at x in the direction dx_n.  This is
                // used both in the predicted reduction as well as the
                // residual predicted reduction. 
                Y_Vector gpxdxn_p_gx;

                // Derivative of the constraint applied to the tangential step
                // this is used in the residual predicted reduction.
                Y_Vector gpxdxt;

                // Norm of gpxdxn_p_gx.  This is used in the penalty parameter
                // computation and predicted reduction. 
                Real norm_gpxdxnpgx;

                // Normal step
                X_Vector dx_n;
                
                // Cauchy point for normal step
                X_Vector dx_ncp;

                // (Corrected) tangential step
                X_Vector dx_t;

                // Tangential step prior to correction
                X_Vector dx_t_uncorrected;
                
                // Cauchy point for tangential step prior to correction
                X_Vector dx_tcp_uncorrected;
                
                // Hessian applied to the normal step.  This is required by
                // W_gradpHdxn as well as the predicted reduction.
                X_Vector H_dxn;

                // Quantity grad f(x) + g'(x)*y + H dx_n projected into the
                // null-space of the constraints.  This is required in the
                // tangential subproblem and the predicted reduction.
                X_Vector W_gradpHdxn;
                
                // Hessian applied to the uncorrected tangential step.  This
                // is needed in the predicted reduction.
                X_Vector H_dxtuncorrected;
                
                // Initialization constructors
                explicit t(X_Vector const & x_,Y_Vector const & y_) : 
                    Unconstrained <Real,XX>::State::t(x_),
                    y(Y::init(y_)),
                    dy(Y::init(y_)),
                    zeta(0.8),
                    eta0(0.5),
                    rho(1.0),
                    rho_old(rho),
                    rho_bar(1e-8),
                    eps_constr(1e-8),
                    xi_qn(1e-4),
                    xi_pg(1e-4),
                    xi_proj(1e-4),
                    xi_tang(1e-4),
                    xi_lmh(1e-4),
                    xi_lmg(1e4),
                    xi_4(2.),
                    rpred(std::numeric_limits<Real>::quiet_NaN()),
                    PSchur_left_type(Operators::Identity),
                    PSchur_right_type(Operators::Identity),
                    augsys_iter_max(100),
                    augsys_rst_freq(0),
                    g_x(Y::init(y_)),
                    norm_gxtyp(std::numeric_limits<Real>::quiet_NaN()),
                    gpxdxn_p_gx(Y::init(y_)),
                    gpxdxt(Y::init(y_)),
                    norm_gpxdxnpgx(std::numeric_limits<Real>::quiet_NaN()),
                    dx_n(X::init(x_)),
                    dx_ncp(X::init(x_)),
                    dx_t(X::init(x_)),
                    dx_t_uncorrected(X::init(x_)),
                    dx_tcp_uncorrected(X::init(x_)),
                    H_dxn(X::init(x_)),
                    W_gradpHdxn(X::init(x_)),
                    H_dxtuncorrected(X::init(x_))
                {
                    Y::copy(y_,y);
                }
            };
            
            // Check that we have a valid set of parameters.  
            static void check_(Messaging const & msg,t const & state) {
                   
                // Use this to build an error message
                std::stringstream ss;
                
                // Check that the fraction of the trust-region used for the
                // quasi-Newton step is between 0 and 1
                // is positive
                if(state.zeta <= Real(0.) || state.zeta >= Real(1.)) 
                    ss << "The fraction of the trust-region used for the "
                        "quasi-Newton step must lie in the interval (0,1): "
                        "zeta = " << state.zeta;
                
                // Check that the trust-region parameter that bounds the
                // error in the preduction reduction lies between 0 and 1-eta1.
                else if(state.eta0 <= Real(0.)
                    || state.eta0 >= Real(1.)-state.eta1) 
                    ss << "The trust-region parameter that bounds the error "
                        "in the predicted reduction must lie in the interval "
                        "(0,1-eta1): eta0 = " << state.eta0;

                // Check that the augmented Lagrangian penalty parameter is
                // greater than or equal to 1
                else if(state.rho < Real(1.))
                    ss << "The augmented Lagrangian penalty parameter must be "
                        "greater than or equal to 1: rho = " << state.rho;

                // Check that the last penalty parameter is greater than or
                // equal to 1
                else if(state.rho_old < Real(1.))
                    ss << "The previous augmented Lagrangian penalty parameter"
                        "must be greater than or equal to 1: rho_old = "
                        << state.rho_old;

                // Check that the fixed increased to the augmented Lagrangian
                // penalty parameter is positive 
                else if(state.rho_bar <= Real(0.))
                    ss << "The fixed increase to the augmented Lagrangian "
                        "penalty paramter must be positive: rho_bar = " 
                        << state.rho_bar;

                // Check that the stopping tolerance for the norm of the
                // constraints is positive
                else if(state.eps_constr <= Real(0.))
                    ss << "The tolerance used in the norm of the constraints "
                        "stopping condition must be positive: eps_constr = "
                        << state.eps_constr;

                // Check that the quasi-Newton step inexactness tolerance lies 
                // in the interval (0,1) 
                else if(state.xi_qn <= Real(0.) || state.xi_qn >= Real(1.))
                    ss << "The quasi-Newton step inexactness tolerance must "
                        "lie in the interval (0,1): xi_qn = " << state.xi_qn;
                
                // Check that the projected gradient inexactness tolerance lies 
                // in the interval (0,1) 
                else if(state.xi_pg <= Real(0.) || state.xi_pg >= Real(1.))
                    ss << "The projected gradient inexactness tolerance must "
                        "lie in the interval (0,1): xi_pg = " << state.xi_pg;
                
                // Check that the nullspace projection inexactness tolerance
                // lies in the interval (0,1) 
                else if(state.xi_proj <= Real(0.) || state.xi_proj >= Real(1.))
                    ss << "The nullspace projection inexactness tolerance must "
                        "lie in the interval (0,1): xi_proj = "<< state.xi_proj;
                
                // Check that the tangential step inexactness tolerance
                // lies in the interval (0,1) 
                else if(state.xi_tang <= Real(0.) || state.xi_tang >= Real(1.))
                    ss << "The tangential step inexactness tolerance must "
                        "lie in the interval (0,1): xi_tang = "<< state.xi_tang;
                
                // Check that the Lagrange multiplier inexactness tolerance
                // lies in the interval (0,1) 
                else if(state.xi_lmh <= Real(0.) || state.xi_lmh >= Real(1.))
                    ss << "The Lagrange multiplier inexactness tolerance must "
                        "lie in the interval (0,1): xi_lmh = " << state.xi_lmh;

                // Check that the absolute tolerance on the residual of the
                // Lagrange multiplier solve is positive
                else if(state.xi_lmg <= Real(0.))
                    ss << "The Lagrange multiplier residual tolerance must "
                        "be positive: xi_lmg = " << state.xi_lmg;

                // Check that the tolerance for the error acceptable in
                // the tangential step is greater than 1.
                else if(state.xi_4 <= Real(1.))
                    ss << "The tolerance on the acceptable error in the "
                        "tangential step must be greater than or equal to 1: "
                        "xi_4 = " << state.xi_4;
                
                // Check that the left preconditioner for the augmented system
                // is either defined by the user or the identity.
                else if(state.PSchur_left_type != Operators::Identity &&
                    state.PSchur_left_type != Operators::UserDefined)
                    ss << "The left preconditioner for the augmented system "
                        "must be either user defined or the identity: "
                        "PSchur_left_type = "
                        << Operators::to_string(state.PSchur_left_type);
                
                // Check that the right preconditioner for the augmented system
                // is either defined by the user or the identity.
                else if(state.PSchur_right_type != Operators::Identity &&
                    state.PSchur_right_type != Operators::UserDefined)
                    ss << "The right preconditioner for the augmented system "
                        "must be either user defined or the identity: "
                        "PSchur_right_type = "
                        << Operators::to_string(state.PSchur_right_type);

                // Check that the number of iterations used when solving the
                // augmented system is positive
                else if(state.augsys_iter_max <= 0)
                    ss << "The number of iterations used when solving the "
                        "augmented system must be positive: augsys_iter_max = "
                        << state.augsys_iter_max;
            }
            static void check(Messaging const & msg,t const & state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {
            // Disallow constructors
            NO_CONSTRUCTORS(Restart);
       
            // Holds restart information
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <Natural> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Y_Vector> > Y_Vectors;

            // Checks whether we have a valid real label.
            struct is_real :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_real()(name) ||
                        name == "zeta" ||
                        name == "eta0" ||
                        name == "rho" ||
                        name == "rho_old" ||
                        name == "rho_bar" ||
                        name == "eps_constr" ||
                        name == "xi_qn" || 
                        name == "xi_pg" ||
                        name == "xi_proj" ||
                        name == "xi_tang" ||
                        name == "xi_lmh" ||
                        name == "xi_lmg" ||
                        name == "xi_4" ||
                        name == "rpred" ||
                        name == "norm_gxtyp" ||
                        name == "norm_gpxdxnpgx" 
                    )
                        return true;
                    else
                        return false;
                    }
            };
            
            // Checks whether we have a valid natural number label.
            struct is_nat :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                        ::is_nat()(name) ||
                        name == "augsys_iter_max" ||
                        name == "augsys_rst_freq"
                    )
                        return true;
                    else
                        return false;
                }
            };
           
            // Checks whether we have a valid parameter label.
            struct is_param :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_param()(name) ||
                        name == "PSchur_right_type" ||
                        name == "PSchur_left_type" 
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid variable label
            struct is_x : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_x()(name) ||
                        name == "dx_n" ||
                        name == "dx_ncp" ||
                        name == "dx_t" ||
                        name == "dx_t_uncorrected" ||
                        name == "dx_tcp_uncorrected" ||
                        name == "H_dxn" ||
                        name == "W_gradpHdxn" ||
                        name == "H_dxtuncorrected" 
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid equality multiplier label
            struct is_y : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( name == "y" ||
                        name == "dy" ||
                        name == "g_x" ||
                        name == "gpxdxn_p_gx" ||
                        name == "gpxdxt"
                    ) 
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have valid labels
            static void checkLabels(
                Messaging const & msg,
                Reals const & reals,
                Nats const & nats,
                Params const & params,
                X_Vectors const & xs,
                Y_Vectors const & ys
            ) {
                Utility::checkLabels(msg,is_real(),reals.first," real name: ");
                Utility::checkLabels(msg,is_nat(),nats.first," natural name: ");
                Utility::checkLabels(
                    msg,is_param(),params.first," paramater name: ");
                Utility::checkLabels(msg,is_x(),xs.first," variable name: ");
                Utility::checkLabels(
                    msg,is_y(),ys.first," equality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string const &,std::string const &,std::string>
            {
                std::string operator() (
                    std::string const & label,
                    std::string const & val
                ) const {

                    // Create a base message
                    const std::string base
                        ="During serialization, found an invalid ";

                    // Used to build the message 
                    std::stringstream ss;

                    // Check the unconstrained parameters
                    if(typename Unconstrained <Real,XX>::Restart
                        ::is_param()(label)
                    ) {
                        ss << typename Unconstrained <Real,XX>::Restart
                            ::checkParamVal()(label,val);

                    // Check the type of the left preconditioner to the
                    // augmented system
                    } else if(label=="PSchur_left_type"){
                        if(!Operators::is_valid()(val))
                            ss << base << "preconditioner type: " << val;
                    
                    // Check the type of the right preconditioner to the
                    // augmented system
                    } else if(label=="PSchur_right_type"){
                        if(!Operators::is_valid()(val))
                            ss << base << "preconditioner type: " << val;
                    }

                    return ss.str();
                }
            };
            
            // Copy out all equality multipliers 
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs,
                Y_Vectors & ys
            ) {
                ys.first.emplace_back("y");
                ys.second.emplace_back(std::move(state.y));
                ys.first.emplace_back("dy");
                ys.second.emplace_back(std::move(state.dy));
                ys.first.emplace_back("g_x");
                ys.second.emplace_back(std::move(state.g_x));
                ys.first.emplace_back("gpxdxn_p_gx");
                ys.second.emplace_back(std::move(state.gpxdxn_p_gx));
                ys.first.emplace_back("gpxdxt");
                ys.second.emplace_back(std::move(state.gpxdxt));
                xs.first.emplace_back("dx_n");
                xs.second.emplace_back(std::move(state.dx_n));
                xs.first.emplace_back("dx_ncp");
                xs.second.emplace_back(std::move(state.dx_ncp));
                xs.first.emplace_back("dx_t");
                xs.second.emplace_back(std::move(state.dx_t));
                xs.first.emplace_back("dx_t_uncorrected");
                xs.second.emplace_back(std::move(state.dx_t_uncorrected));
                xs.first.emplace_back("dx_tcp_uncorrected");
                xs.second.emplace_back(std::move(state.dx_tcp_uncorrected));
                xs.first.emplace_back("H_dxn");
                xs.second.emplace_back(std::move(state.H_dxn));
                xs.first.emplace_back("W_gradpHdxn");
                xs.second.emplace_back(std::move(state.W_gradpHdxn));
                xs.first.emplace_back("H_dxtuncorrected");
                xs.second.emplace_back(std::move(state.H_dxtuncorrected));
            }

            // Copy out all the scalar information
            static void stateToScalars(
                typename State::t & state,
                Reals & reals,
                Nats & nats,
                Params & params
            ) { 
                // Copy in all the real numbers
                reals.first.emplace_back("zeta");
                reals.second.emplace_back(state.zeta);
                reals.first.emplace_back("eta0");
                reals.second.emplace_back(state.eta0);
                reals.first.emplace_back("rho");
                reals.second.emplace_back(state.rho);
                reals.first.emplace_back("rho_old");
                reals.second.emplace_back(state.rho_old);
                reals.first.emplace_back("rho_bar");
                reals.second.emplace_back(state.rho_bar);
                reals.first.emplace_back("eps_constr");
                reals.second.emplace_back(state.eps_constr);
                reals.first.emplace_back("xi_qn");
                reals.second.emplace_back(state.xi_qn);
                reals.first.emplace_back("xi_pg");
                reals.second.emplace_back(state.xi_pg);
                reals.first.emplace_back("xi_proj");
                reals.second.emplace_back(state.xi_proj);
                reals.first.emplace_back("xi_tang");
                reals.second.emplace_back(state.xi_tang);
                reals.first.emplace_back("xi_lmh");
                reals.second.emplace_back(state.xi_lmh);
                reals.first.emplace_back("xi_lmg");
                reals.second.emplace_back(state.xi_lmg);
                reals.first.emplace_back("xi_4");
                reals.second.emplace_back(state.xi_4);
                reals.first.emplace_back("rpred");
                reals.second.emplace_back(state.rpred);
                reals.first.emplace_back("norm_gxtyp");
                reals.second.emplace_back(state.norm_gxtyp);
                reals.first.emplace_back("norm_gpxdxnpgx");
                reals.second.emplace_back(state.norm_gpxdxnpgx);

                // Copy in all the natural numbers
                nats.first.emplace_back("augsys_iter_max");
                nats.second.emplace_back(state.augsys_iter_max);
                nats.first.emplace_back("augsys_rst_freq");
                nats.second.emplace_back(state.augsys_rst_freq);

                // Copy in all the parameters
                params.first.emplace_back("PSchur_left_type");
                params.second.emplace_back(
                    Operators::to_string(state.PSchur_left_type));
                params.first.emplace_back("PSchur_right_type");
                params.second.emplace_back(
                    Operators::to_string(state.PSchur_right_type));
            }
            
            // Copy in all equality multipliers 
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys
            ) {
                typename std::list <X_Vector>::iterator y
                    =ys.second.begin();
                for(typename std::list <std::string>::iterator name 
                        =ys.first.begin();
                    name!=ys.first.end();
                ) {
                    // Make a copy of the current iterators.  We use these
                    // to remove elements
                    typename std::list <std::string>::iterator name0 = name;
                    typename std::list <Y_Vector>::iterator y0 = y;

                    // Increment our primary iterators 
                    name++; y++;

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name0=="y")
                        state.y = std::move(*y0);
                    else if(*name0=="dy")
                        state.dy = std::move(*y0);
                    else if(*name0=="g_x")
                        state.g_x = std::move(*y0);
                    else if(*name0=="gpxdxn_p_gx")
                        state.gpxdxn_p_gx = std::move(*y0);
                    else if(*name0=="gpxdxt")
                        state.gpxdxt = std::move(*y0);
                }

                typename std::list <X_Vector>::iterator x
                    =xs.second.begin();
                for(typename std::list <std::string>::iterator name 
                        =xs.first.begin();
                    name!=xs.first.end();
                ) {
                    // Make a copy of the current iterators.  We use these
                    // to remove elements
                    typename std::list <std::string>::iterator name0 = name;
                    typename std::list <X_Vector>::iterator x0 = x;

                    // Increment our primary iterators 
                    name++; x++;

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name0=="dx_n")
                        state.dx_n = std::move(*x0);
                    else if(*name0=="dx_ncp")
                        state.dx_ncp = std::move(*x0);
                    else if(*name0=="dx_t")
                        state.dx_t = std::move(*x0);
                    else if(*name0=="dx_t_uncorrected")
                        state.dx_t_uncorrected = std::move(*x0);
                    else if(*name0=="dx_tcp_uncorrected")
                        state.dx_tcp_uncorrected = std::move(*x0);
                    else if(*name0=="H_dxn")
                        state.H_dxn = std::move(*x0);
                    else if(*name0=="W_gradpHdxn")
                        state.W_gradpHdxn = std::move(*x0);
                    else if(*name0=="H_dxtuncorrected")
                        state.H_dxtuncorrected = std::move(*x0);
                }
            }
            
            // Copy in all the scalar information
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Nats & nats,
                Params & params
            ) { 
                // Copy in any reals 
                typename std::list <Real>::iterator real=reals.second.begin();
                for(std::list <std::string>::iterator name=reals.first.begin();
                    name!=reals.first.end();
                    name++,real++
                ){
                    if(*name=="zeta") state.zeta=*real;
                    else if(*name=="eta0") state.eta0=*real;
                    else if(*name=="rho") state.rho=*real;
                    else if(*name=="rho_old") state.rho_old=*real;
                    else if(*name=="rho_bar") state.rho_bar=*real;
                    else if(*name=="eps_constr") state.eps_constr=*real;
                    else if(*name=="xi_qn") state.xi_qn=*real;
                    else if(*name=="xi_pg") state.xi_pg=*real;
                    else if(*name=="xi_proj") state.xi_proj=*real;
                    else if(*name=="xi_tang") state.xi_tang=*real;
                    else if(*name=="xi_lmh") state.xi_lmh=*real;
                    else if(*name=="xi_lmg") state.xi_lmg=*real;
                    else if(*name=="xi_4") state.xi_4=*real;
                    else if(*name=="rpred") state.rpred=*real;
                    else if(*name=="norm_gxtyp") state.norm_gxtyp=*real;
                    else if(*name=="norm_gpxdxnpgx") state.norm_gpxdxnpgx=*real;
                }
                
                // Next, copy in any naturals
                std::list <Natural>::iterator nat=nats.second.begin();
                for(std::list <std::string>::iterator name=nats.first.begin();
                    name!=nats.first.end();
                    name++,nat++
                ){
                    if(*name=="augsys_iter_max") state.augsys_iter_max=*nat;
                    else if(*name=="augsys_rst_freq")state.augsys_rst_freq=*nat;
                }
                
                // Next, copy in any parameters 
                std::list <std::string>::iterator param=params.second.begin();
                for(std::list <std::string>::iterator name=params.first.begin();
                    name!=params.first.end();
                    name++,param++
                ){
                    if(*name=="PSchur_left_type")
                        state.PSchur_left_type=Operators::from_string(*param);
                    else if(*name=="PSchur_right_type")
                        state.PSchur_right_type=Operators::from_string(*param);
                }
            }

            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>
                    ::Restart::stateToVectors(state,xs);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::stateToVectors(state,xs,ys);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::stateToScalars(state,reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::stateToScalars(state,reals,nats,params);
            }

            // Capture data from structures controlled by the user.  
            static void capture(
                Messaging const & msg,
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,ys);

                // Check the strings used to represent parameters
                Utility::checkParams(msg,checkParamVal(),params);

                // Copy in the variables 
                Unconstrained <Real,XX>
                    ::Restart::vectorsToState(state,xs);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::vectorsToState(state,xs,ys);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::scalarsToState(state,reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::scalarsToState(state,reals,nats,params);

                // Check that we have a valid state 
                State::check(msg,state);
            }
        };
        
        // All the functions required by an optimization algorithm.  Note, this
        // routine owns the memory for these operations.  
        struct Functions {
            // Disallow constructors
            NO_CONSTRUCTORS(Functions);
            
            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t);

                // Equality constraints 
                std::unique_ptr <VectorValuedFunction <Real,XX,YY> > g;

                // Left preconditioner for the augmented system
                std::unique_ptr <Operator <Real,YY,YY> > PSchur_left;

                // Right preconditioner for the augmented system
                std::unique_ptr <Operator <Real,YY,YY> > PSchur_right;
                
                // Initialize all of the pointers to null
                t() : Unconstrained <Real,XX>::Functions::t(), g(nullptr),
                    PSchur_left(nullptr), PSchur_right(nullptr) {}
            };

            struct EqualityModifications
                : public Optizelle::ScalarValuedFunctionModifications <Real,XX>
            {
            public:
                // Disallow constructors
                NO_COPY_ASSIGNMENT(EqualityModifications); 

            private:
                // Underlying modification.  This takes control of the memory
                std::unique_ptr <
                    Optizelle::ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Equality constraint.
                Optizelle::VectorValuedFunction <Real,XX,YY> const & g;

                // Reference to equality Lagrange multiplier
                Y_Vector const & y;

                // Reference to parameter for the augmented-Lagrangian
                Real const & rho;

                // Some workspace for the below functions
                mutable X_Vector grad_tmp;
                mutable X_Vector x_tmp1;
                mutable Y_Vector y_tmp1;

                // Variables used for caching.  The boolean values denote
                // whether or not we've started caching yet.
                mutable std::pair <bool,X_Vector> x_merit;
                mutable Y_Vector g_x;
                mutable std::pair <bool,X_Vector> x_grad;
                mutable std::pair <bool,Y_Vector> y_grad;
                mutable X_Vector gpxsy; 

                // Adds the Lagrangian pieces to the gradient
                void grad_lag(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_lag
                ) const {
                    // grad_lag <- grad f(x)
                    X::copy(grad,grad_lag);
                    
                    // If relative error between the current and cached values
                    // is large, compute anew.
                    if( rel_err_cached <Real,XX> (x,x_grad)
                            >= std::numeric_limits <Real>::epsilon()*1e1 ||
                        rel_err_cached <Real,YY> (y,y_grad)
                            >= std::numeric_limits <Real>::epsilon()*1e1 
                    ) {
                        // gpxsy <- g'(x)* y 
                        g.ps(x,y,gpxsy);

                        // Cache the values
                        x_grad.first=true;
                        X::copy(x,x_grad.second);
                        y_grad.first=true;
                        Y::copy(y,y_grad.second);
                    }

                    // grad <- grad f(x) + g'(x)*y 
                    X::axpy(Real(1.),gpxsy,grad_lag);
                }

            public:
                EqualityModifications(
                    typename State::t const & state,
                    typename Functions::t & fns
                ) : f_mod(std::move(fns.f_mod)),
                    g(*(fns.g)),
                    y(state.y),
                    rho(state.rho),
                    grad_tmp(X::init(state.x)),
                    x_tmp1(X::init(state.x)),
                    y_tmp1(Y::init(state.y)),
                    x_merit(false,X::init(state.x)),
                    g_x(Y::init(state.y)),
                    x_grad(false,X::init(state.x)),
                    y_grad(false,Y::init(state.y)),
                    gpxsy(X::init(state.x))
                { }

                // Merit function additions to the objective
                virtual Real merit(X_Vector const & x,Real const & f_x) const {
                    // Do the underlying modification of the objective
                    Real merit_x = f_mod->merit(x,f_x);
                    
                    // If we've not started caching or the relative error
                    // is large, compute anew.
                    if( rel_err_cached <Real,XX> (x,x_merit)
                            >= std::numeric_limits <Real>::epsilon()*1e1
                    ) {
                        // g_x <- g(x)
                        g(x,g_x);
                    
                        // Cache the values
                        x_merit.first=true;
                        X::copy(x,x_merit.second);
                    }

                    // Return f(x) + < y,g(x) > + rho || g(x) ||^2   
                    return merit_x + Y::innr(y,g_x) + rho * Y::innr(g_x,g_x);
                }

                // Stopping condition modification of the gradient
                virtual void grad_stop(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_stop
                ) const {
                    f_mod->grad_stop(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_stop);
                }

                // Diagnostic modification of the gradient
                virtual void grad_diag(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_diag
                ) const {
                    f_mod->grad_diag(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_diag);
                }

                // Modification of the gradient when finding a trial step
                virtual void grad_step(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_step
                ) const {
                    f_mod->grad_step(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_step);
                }

                // Modification of the gradient for a quasi-Newton method 
                virtual void grad_quasi(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_quasi
                ) const {
                    f_mod->grad_quasi(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_quasi);
                }

                // Modification of the gradient when solving for the equality
                // multiplier
                virtual void grad_mult(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_mult
                ) const {
                    f_mod->grad_mult(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_mult);
                }

                // Modification of the Hessian-vector product when finding a
                // trial step
                virtual void hessvec_step(
                    X_Vector const & x,
                    X_Vector const & dx,
                    X_Vector const & H_dx,
                    X_Vector & Hdx_step 
                ) const {
                    // Modify the Hessian vector product 
                    f_mod->hessvec_step(x,dx,H_dx,Hdx_step);

                    // x_tmp1 <- (g''(x)dx)*y
                    g.pps(x,dx,y,x_tmp1);

                    // Hdx_step <- hess f(x)dx + (g''(x)dx)*y  
                    X::axpy(Real(1.),x_tmp1,Hdx_step);
                }
            };

            // The identity operator 
            struct Identity : public Operator <Real,YY,YY> {
                void operator () (Y_Vector const & dy,Y_Vector & result) const{
                    Y::copy(dy,result);
                }
            };

            // Check that all the functions are defined
            static void check(Messaging const & msg,t const & fns) {

                // Check the unconstrained pieces
                Unconstrained <Real,XX>::Functions::check(msg,fns);
                
                // Check that the equality constraints exist 
                if(fns.g.get()==nullptr)
                    msg.error("Missing the equality constraint definition.");

                // Check that preconditioners exist
                if(fns.PSchur_left.get()==nullptr)
                    msg.error("Missing a left preconditioner for the "
                        "augmented system.");
                if(fns.PSchur_right.get()==nullptr)
                    msg.error("Missing a right preconditioner for the "
                        "augmented system.");
            }

            // Initialize any missing functions for just equality constrained 
            // optimization.
            static void init_(
                Messaging const & msg,
                typename State::t const & state,
                t & fns
            ) {
                // Determine the left preconditioner for the augmented system
                switch(state.PSchur_left_type){
                    case Operators::Identity:
                        fns.PSchur_left.reset(new Identity());
                        break;
                    case Operators::UserDefined:
                        if(fns.PSchur_left.get()==nullptr)
                            msg.error("An externally defined left "
                                "preconditioner for the augmented system must "
                                "be provided explicitly.");
                        break;
                    default:
                        msg.error("Not a valid left preconditioner for the "
                            "augmented system.");
                        break;
                }
                
                // Determine the right preconditioner for the augmented system
                switch(state.PSchur_right_type){
                    case Operators::Identity:
                        fns.PSchur_right.reset(new Identity());
                        break;
                    case Operators::UserDefined:
                        if(fns.PSchur_right.get()==nullptr)
                            msg.error("An externally defined right "
                                "preconditioner for the augmented system must "
                                "be provided explicitly.");
                        break;
                    default:
                        msg.error("Not a valid right preconditioner for the "
                            "augmented system.");
                        break;
                }

                // If a trust-region operator has not been provided, use the
                // identity.
                if(fns.TRS.get()==nullptr)
                    fns.TRS.reset(new typename
                        Unconstrained <Real,XX>::Functions::Identity());
                
                // Check that all functions are defined 
                check(msg,fns);
                
                // Modify the objective 
                fns.f_mod.reset(new EqualityModifications(state,fns));
            }

            // Initialize any missing functions 
            static void init(
                Messaging const & msg,
                typename State::t const & state,
                t & fns
            ) {
                Unconstrained <Real,XX>
                    ::Functions::init_(msg,state,fns);
                EqualityConstrained <Real,XX,YY>
                    ::Functions::init_(msg,state,fns);
            }
        };
        
        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
            // Disallow constructors
            NO_CONSTRUCTORS(Printer);

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { 
                // Norm of the constrained 
                out.emplace_back(Utility::atos("||g(x)||"));
                    
                // Trust-region information
                out.emplace_back(Utility::atos("ared"));
                out.emplace_back(Utility::atos("pred"));
                out.emplace_back(Utility::atos("ared/pred"));
                   
                // Krylov method information
                out.emplace_back(Utility::atos("KryIter"));
                out.emplace_back(Utility::atos("KryErr"));
                out.emplace_back(Utility::atos("KryWhy"));
            }
            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer::getStateHeader_
                    (state,out);
            }

            // Gets the state information for output
            static void getState_(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                // Create some shortcuts
                Y_Vector const & g_x = state.g_x;
                Natural const & krylov_iter=state.krylov_iter;
                Real const & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Natural const & iter=state.iter;
                Natural const & rejected_trustregion=state.rejected_trustregion;
                Real const & pred = state.pred;
                Real const & ared = state.ared;

                // Figure out if we're at the absolute beginning of the
                // optimization.  We have to be a little saavy about this
                // since we could be on the first iteration, but in the
                // middle of rejecting trust-region steps and still want 
                // to output things.
                bool opt_begin = (iter==1) &&
                        (rejected_trustregion == 0);

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Norm of the gradient 
                Real norm_gx = sqrt(Y::innr(g_x,g_x));
                out.emplace_back(Utility::atos(norm_gx));
                    
                // Actual vs. predicted reduction 
                if(!opt_begin) {
                    out.emplace_back(Utility::atos(ared));
                    out.emplace_back(Utility::atos(pred));
                    out.emplace_back(Utility::atos(ared/pred));
                } else 
                    for(Natural i=0;i<3;i++)
                        out.emplace_back(Utility::blankSeparator);
                
                // Krylov method information
                if(!opt_begin) {
                    out.emplace_back(Utility::atos(krylov_iter));
                    out.emplace_back(Utility::atos(krylov_rel_err));
                    out.emplace_back(Utility::atos(krylov_stop));
                } else 
                    for(Natural i=0;i<3;i++)
                        out.emplace_back(Utility::blankSeparator);

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=Utility::blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                bool const & noiter,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getState_(fns,state,blank,out);
            }
            
            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) { }

            // Combines all of the Krylov information
            static void getKrylov(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };

        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms);

            // The operator for the augmented system,
            //
            // [ I      g'(x)* ]
            // [ g'(x)  0      ]
            //
            struct AugmentedSystem: public Operator <Real,XXxYY,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
                X_Vector const & x_base;
            public:
                AugmentedSystem(
                    typename State::t const & state_,
                    typename Functions::t const & fns_,
                    X_Vector const & x_base_
                ) : state(state_), fns(fns_), x_base(x_base_) {}
                
                // Operator interface
                void operator () (
                    const XxY_Vector & dx_dy,
                    XxY_Vector & result
                ) const{
                    // Create some shortcuts
                    VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);

                    // g'(x_base)* dy
                    g.ps(x_base,dx_dy.second,result.first);

                    // dx + g'(x_base)* dy 
                    X::axpy(Real(1.),dx_dy.first,result.first);

                    // g'(x_base)* dx
                    g.p(x_base,dx_dy.first,result.second);
                }
            };
            
            // The block diagonal preconditioner 
            //
            // [ PH_x    0    ]
            // [ 0       PH_y ]
            //
            struct BlockDiagonalPreconditioner:
                public Operator <Real,XXxYY,XXxYY> {
            private:
                Operator <Real,XX,XX> const & PH_x;
                const Operator <Real,YY,YY>& PH_y;
            public:
                BlockDiagonalPreconditioner(
                    Operator <Real,XX,XX> const & PH_x_,
                    const Operator <Real,YY,YY>& PH_y_ 
                ) : PH_x(PH_x_), PH_y(PH_y_) {}
                
                // Operator interface
                void operator () (
                    const XxY_Vector & dx_dy,
                    XxY_Vector & result
                ) const{
                    // PH_x dx
                    PH_x(dx_dy.first,result.first);
                    
                    // PH_y dy
                    PH_y(dx_dy.second,result.second);
                }
            };

            // Sets the tolerances for the quasi-normal Newton solve
            struct QNManipulator : GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                explicit QNManipulator(
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_qn = state.xi_qn;
                    Real const & norm_gxtyp = state.norm_gxtyp;
                    Real const & eps_constr= state.eps_constr;

                    // Find || g'(x)dx_ncp + g(x) ||
                    Real norm_gpxdxncp_p_g = sqrt(Y::innr(bb.second,bb.second));

                    // Return xi_qn * || g'(x)dx_ncp + g(x) ||
                    eps = xi_qn * norm_gpxdxncp_p_g;

                    // If the Cauchy point actually brings us to optimality,
                    // it's hard to hit the tolerance above.  In this case,
                    // try to detect the condition and bail early.  The way
                    // we detect this is by noting that
                    //
                    // g(x+dx) ~= g(x)+g'(x)dx
                    //
                    // Hence, if || g(x)+g'(x)dx || is small, we should be
                    // feasible after the step.  Therefore, we check if we
                    // satisfy our stopping condition for feasibility.  If so,
                    // we bail.
                    if(norm_gpxdxncp_p_g < eps_constr*norm_gxtyp)
                        eps=Real(1.);
                }
            };

            // Finds the quasi-normal step
            static void quasinormalStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);
                X_Vector const & x=state.x;
                Y_Vector const & g_x=state.g_x;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                Real const & delta = state.delta;
                Real const & zeta = state.zeta;
                X_Vector & dx_ncp=state.dx_ncp;
                X_Vector & dx_n=state.dx_n;

                // Find the Cauchy point.

                // Find g'(x)*g(x)
                X_Vector gps_g(X::init(x));
                g.ps(x,g_x,gps_g);

                // Find g'(x)g'(x)*g(x)
                Y_Vector gp_gps_g(Y::init(g_x));
                g.p(x,gps_g,gp_gps_g);

                // Find || g'(x)*g(x) ||^2
                Real norm_gpsg_2 = X::innr(gps_g,gps_g);

                // Find || g'(x)g'(x)*g(x) ||^2
                Real norm_gpgpsg_2 = Y::innr(gp_gps_g,gp_gps_g);

                // Find the Cauchy point,
                // -|| g'(x)*g(x) ||^2 / || g'(x)g'(x)*g(x) ||^2 g'(x)*g(x)
                X::copy(gps_g,dx_ncp);
                X::scal(-norm_gpsg_2/norm_gpgpsg_2,dx_ncp);

                // If || dx_ncp || >= zeta delta, scale it back to zeta
                // delta and return
                Real norm_dxncp = sqrt(X::innr(dx_ncp,dx_ncp));
                if(norm_dxncp >= zeta*delta) {
                    X::scal(zeta*delta/norm_dxncp,dx_ncp);
                    X::copy(dx_ncp,dx_n);
                    return;
                }

                // Find the Newton step

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(g_x));
                XxY::zero(x0);

                // Create the rhs, b0=(-dx_ncp,-g'(x)dx_ncp-g(x)) 
                XxY_Vector b0(XxY::init(x0));
                X::copy(dx_ncp,b0.first);
                X::scal(Real(-1.),b0.first);
                g.p(x,dx_ncp,b0.second);
                Y::scal(Real(-1.),b0.second);
                Y::axpy(Real(-1.),g_x,b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l (I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r (I,*(fns.PSchur_right));

                // Solve the augmented system for the Newton step
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    QNManipulator(state,fns),
                    x0 
                );

                // Find the Newton shift, dx_dnewton = dx_newton-dx_ncp
                X_Vector & dx_dnewton = x0.first;

                // Find the Newton step
                X::copy(dx_ncp,dx_n);
                X::axpy(Real(1.),dx_dnewton,dx_n);

                // If the dx_n is smaller than zeta deta, then return
                // it as the quasi-normal step
                Real norm_dxnewton = sqrt(X::innr(dx_n,dx_n));
                if(norm_dxnewton <= zeta*delta) return;

                // Otherwise, compute the dogleg step.  In order to accomplish
                // this, we need to find theta so that
                // || dx_ncp + theta dx_dnewton || = zeta*delta
                // and then set dx_n = dx_ncp + theta dx_dnewton.
                Real aa = X::innr(dx_dnewton,dx_dnewton);
                Real bb = Real(2.) * X::innr(dx_dnewton,dx_ncp);
                Real cc = norm_dxncp*norm_dxncp - zeta*zeta*delta*delta;
                Natural nroots;
                Real r1;
                Real r2;
                quad_equation(aa,bb,cc,nroots,r1,r2);
                Real theta = r1 > r2 ? r1 : r2;
                X::copy(dx_ncp,dx_n);
                X::axpy(theta,dx_dnewton,dx_n);
            }
            
            // Sets the tolerances for projecting 
            //
            // grad f(x) + g'(x)*y + H dx_n
            //
            // into the null space of g'(x).
            struct NullspaceProjForGradLagPlusHdxnManipulator
                : GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                explicit NullspaceProjForGradLagPlusHdxnManipulator(
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_pg = state.xi_pg;
                    Real const & delta = state.delta;

                    // Find || W (grad L(x,y) + H dx_n) || = || xx_1 || 
                    Real norm_WgpHdxn = sqrt(X::innr(xx.first,xx.first));

                    // Find || grad L(x,y) + H dx_n || = || bb_1 ||
                    Real norm_gradpHdxn = sqrt(X::innr(bb.first,bb.first));

                    // The bound is xi_pg min( || W (grad L(x,y) + H dx_n) ||,
                    // delta, || grad L(x,y) + H dx_n || )
                    eps = norm_WgpHdxn < delta ? norm_WgpHdxn : delta;
                    eps = eps < norm_gradpHdxn ? eps : norm_gradpHdxn;
                    eps = xi_pg*eps;

                    // If the projected gradient is in the nullspace of
                    // the constraints, it's hard to hit the tolerance above.
                    // In this case, try to detect the condition and bail early.
                    // Specifically, sometimes the right hand side is
                    // a decent size, but the solution is small.  Therefore,
                    // we exit when the norm of the current projected gradient
                    // is smaller than the norm of the unprojected gradient
                    // by two orders of magnitude larger than machine precision.
                    if(iter >= 2
                        && norm_WgpHdxn 
                            < std::numeric_limits <Real>::epsilon()
                              * norm_gradpHdxn * Real(1e2)
                    )
                        eps=Real(1.);
                }
            };

            // Projects the quantity:
            //
            // grad f(x) + g'(x)*y + H dx_n
            //
            // into the null space of g'(x).  This is required for the
            // tangential subproblem as well as the predicted reduction.
            // Note, this also computes and caches H dx_n.
            static void projectedGradLagrangianPlusHdxn(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & grad=state.grad;
                X_Vector const & H_dxn=state.H_dxn;
                Y_Vector const & y=state.y;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                X_Vector & W_gradpHdxn=state.W_gradpHdxn;

                // Find the gradient modifications for the step computation
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
               
                // Add the Hessian modifications to H(x)dx_n
                X_Vector Hdxn_step(X::init(x));
                    f_mod.hessvec_step(x,dx_n,H_dxn,Hdxn_step);

                // grad_p_Hdxn <- H dxn_step
                X_Vector grad_p_Hdxn(X::init(x));
                    X::copy(Hdxn_step,grad_p_Hdxn);

                // grad_p_Hdxn <- grad f(x) + H dx_n
                X::axpy(Real(1.),grad_step,grad_p_Hdxn);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(y));
                    XxY::zero(x0);

                // Create the rhs, b0=(grad f(x) + H dx_n,0)
                XxY_Vector b0(XxY::init(x0));
                    X::copy(grad_p_Hdxn,b0.first);
                    Y::zero(b0.second);
            
                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l (I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r (I,*(fns.PSchur_right));

                // Solve the augmented system for the nullspace projection 
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    NullspaceProjForGradLagPlusHdxnManipulator(state,fns),
                    x0 
                );

                // Copy out the solution
                X::copy(x0.first,W_gradpHdxn);
            }
            
            // Sets the tolerances for the nullspace projector that projects
            // the current direction in the projected Krylov method. 
            struct NullspaceProjForKrylovMethodManipulator
                : GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                explicit NullspaceProjForKrylovMethodManipulator (
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_proj = state.xi_proj;

                    // Find || W dx_t_uncorrected || = || xx_1 || 
                    Real norm_Wdxt_uncorrected
                        = sqrt(X::innr(xx.first,xx.first));

                    // Find || dx_t_uncorrected || = || bb_1 ||
                    Real norm_dxt_uncorrected
                        = sqrt(X::innr(bb.first,bb.first));

                    // The bound is xi_proj min( || W dx_t_uncorrected  ||,
                    // || dx_t_uncorrected || )
                    eps = norm_Wdxt_uncorrected < norm_dxt_uncorrected
                        ? norm_Wdxt_uncorrected : norm_dxt_uncorrected;
                    eps = xi_proj*eps; 

                    // If the projected direction is in the nullspace of
                    // the constraints, it's hard to hit the tolerance above.
                    // In this case, try to detect the condition and bail early.
                    // Specifically, sometimes the right hand side is
                    // a decent size, but the solution is small.  Therefore,
                    // we exit when the norm of the projected Krylov iterate 
                    // is smaller than the norm of the unprojected iterate 
                    // by two orders of magnitude larger than machine precision.
                    if(iter >= 2
                        && norm_Wdxt_uncorrected
                            < std::numeric_limits <Real>::epsilon()
                              * norm_dxt_uncorrected * Real(1e2)
                    )
                        eps=Real(1.);
                }
            };
            
            // Nullspace projector that projects the current direction in the
            // projected Krylov method. 
            struct NullspaceProjForKrylovMethod: public Operator <Real,XX,XX> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                NullspaceProjForKrylovMethod(
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
               
                // Project dx_t into the nullspace of g'(x)
                void operator () (
                    X_Vector const & dx_t_uncorrected,
                    X_Vector & result
                ) const{
                    // Create some shortcuts
                    X_Vector const & x=state.x;
                    Y_Vector const & y=state.y;
                    unsigned int const & augsys_iter_max=state.augsys_iter_max;
                    unsigned int const & augsys_rst_freq=state.augsys_rst_freq;

                    // Create the initial guess, x0=(0,0)
                    XxY_Vector x0(X::init(x),Y::init(y));
                        XxY::zero(x0);

                    // Create the rhs, b0=(dx_t_uncorrected,0)
                    XxY_Vector b0(XxY::init(x0));
                        X::copy(dx_t_uncorrected,b0.first);
                        Y::zero(b0.second);
                
                    // Build Schur style preconditioners
                    typename Unconstrained <Real,XX>::Functions::Identity I;
                    BlockDiagonalPreconditioner
                        PAugSys_l(I,*(fns.PSchur_left));
                    BlockDiagonalPreconditioner
                        PAugSys_r(I,*(fns.PSchur_right));

                    // Solve the augmented system for the nullspace projection 
                    Optizelle::gmres <Real,XXxYY> (
                        AugmentedSystem(state,fns,x),
                        b0,
                        Real(1.), // This will be overwritten by the manipulator
                        augsys_iter_max,
                        augsys_rst_freq,
                        PAugSys_l,
                        PAugSys_r,
                        NullspaceProjForKrylovMethodManipulator(state,fns),
                        x0 
                    );

                    // Copy out the solution
                    X::copy(x0.first,result);
                }
            };
            
            // Solves the tangential subproblem 
            static void tangentialSubProblem(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & W_gradpHdxn=state.W_gradpHdxn;
                Real const & delta = state.delta;
                Real const & eps_krylov=state.eps_krylov;
                Natural const & krylov_iter_max=state.krylov_iter_max;
                Natural const & krylov_orthog_max=state.krylov_orthog_max;
                KrylovSolverTruncated::t const & krylov_solver
                    = state.krylov_solver;
                X_Vector & dx_t_uncorrected=state.dx_t_uncorrected;
                X_Vector & dx_tcp_uncorrected=state.dx_tcp_uncorrected;
                Real & krylov_rel_err=state.krylov_rel_err;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                
                // Create shortcuts to the functions that we need
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                Operator <Real,XX,XX> const & TRS=*(fns.TRS);
                    
                // Setup the Hessian operator and allocate memory for the
                // Cauchy point.
                typename Unconstrained <Real,XX>::Algorithms::HessianOperator
                    H(f,f_mod,x);

                // Find the quantity - W (g + H dxn).  We use this as the
                // RHS in the linear system solve.
                X_Vector minus_W_gradpHdxn(X::init(x));
                    X::copy(W_gradpHdxn,minus_W_gradpHdxn);
                    X::scal(Real(-1.),minus_W_gradpHdxn);

                // Keep track of the residual errors
                Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                Real residual_err(std::numeric_limits <Real>::quiet_NaN());
            
                switch(krylov_solver) {
                // Truncated conjugate direction
                case KrylovSolverTruncated::ConjugateDirection:
                    truncated_cd(
                        H,
                        minus_W_gradpHdxn,
                        NullspaceProjForKrylovMethod(state,fns), // Add in PH?
                        TRS,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        delta,
                        dx_n,
			true,
                        dx_t_uncorrected,
                        dx_tcp_uncorrected,
                        residual_err0,
                        residual_err,
                        krylov_iter,
                        krylov_stop);
                    break;

                // Truncated MINRES 
                case KrylovSolverTruncated::MINRES:
                    truncated_minres(
                        H,
                        minus_W_gradpHdxn,
                        NullspaceProjForKrylovMethod(state,fns), // Add in PH?
                        TRS,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        delta,
                        dx_n,
                        dx_t_uncorrected,
                        dx_tcp_uncorrected,
                        residual_err0,
                        residual_err,
                        krylov_iter,
                        krylov_stop);

                    // Force a descent direction
                    if(X::innr(dx_t_uncorrected,W_gradpHdxn) > 0)
                        X::scal(Real(-1.),dx_t_uncorrected);
                    if(X::innr(dx_tcp_uncorrected,W_gradpHdxn) > 0)
                        X::scal(Real(-1.),dx_tcp_uncorrected);
                    break;
                }
                krylov_rel_err = residual_err 
                    / (std::numeric_limits <Real>::epsilon()+residual_err0);
                krylov_iter_total += krylov_iter;
            }
            
            // Sets the tolerances for the computation of the tangential
            // step.
            struct TangentialStepManipulator
                : GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                TangentialStepManipulator (
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    X_Vector const & dx_n=state.dx_n;
                    Real const & xi_tang = state.xi_tang;
                    Real const & delta = state.delta;

                    // dxn_p_dxt <- dx_n + dx_t
                    X_Vector dxn_p_dxt(X::init(dx_n));
                    X::copy(dx_n,dxn_p_dxt);
                    X::axpy(Real(1.),xx.first,dxn_p_dxt);

                    // Find || dx_n + dx_t || 
                    Real norm_dxnpdxt = sqrt(X::innr(dxn_p_dxt,dxn_p_dxt));

                    // Find || dx_t_uncorrected || = || bb_1 || 
                    Real norm_dxt_uncorrected= sqrt(X::innr(bb.first,bb.first));

                    // The bound is
                    // delta * min( delta, || dx_n + dx_t ||,
                    //      xi_tang ||dx_t_uncorrected||/delta)
                    eps = delta < norm_dxnpdxt ?  delta : norm_dxnpdxt;
                    eps = eps < xi_tang*norm_dxt_uncorrected/delta
                        ? eps : xi_tang*norm_dxt_uncorrected/delta;
                    eps = eps*delta;
                }
            };
            
            // Finds the tangential step 
            static void tangentialStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                X_Vector const & x=state.x;
                Y_Vector const & y=state.y;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                X_Vector const & dx_t_uncorrected=state.dx_t_uncorrected;
                X_Vector & dx_t=state.dx_t;

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(y));
                    XxY::zero(x0);

                // Create the rhs, b0=(dx_t_uncorrected,0);
                XxY_Vector b0(XxY::init(x0));
                    X::copy(dx_t_uncorrected,b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the tangential step 
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    TangentialStepManipulator(state,fns),
                    x0 
                );

                // Copy out the tangential step
                X::copy(x0.first,dx_t);
            }
            
            // Sets the tolerances for the computation of the Lagrange 
            // multiplier.
            struct LagrangeMultiplierStepManipulator
                : GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t const & state;
                typename Functions::t const & fns;
            public:
                LagrangeMultiplierStepManipulator (
                    typename State::t const & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_lmh = state.xi_lmh;
                    Real const & xi_lmg = state.xi_lmg;
                
                    // Find the norm of the gradient of the Lagrangian.
                    // Sometimes, this is -grad L(x+dx,y).  Sometimes, this
                    // is -grad L(x,y).  In both cases, we just look at the
                    // first element of the RHS.
                    Real norm_grad = sqrt(X::innr(bb.first,bb.first));

                    // The bound is
                    // min( xi_lmg, xi_lmh || grad f(x) + g'(x)*y ||)
                    eps = xi_lmg < norm_grad*xi_lmh ? xi_lmg : norm_grad*xi_lmh;
                }
            };

            // Finds the Lagrange multiplier at the current iterate 
            static void lagrangeMultiplier(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                X_Vector const & grad=state.grad;
                Y_Vector & y=state.y;

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_mult(X::init(grad));
                    f_mod.grad_mult(x,grad,grad_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(y));
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x,y),0);
                XxY_Vector b0(XxY::init(x0));
                    X::copy(grad_mult,b0.first);
                    X::scal(Real(-1.),b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the initial Lagrange
                // multiplier 
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    LagrangeMultiplierStepManipulator(state,fns),
                    x0 
                );

                // Find the Lagrange multiplier based on this step
                Y::axpy(Real(1.),x0.second,y);
            }
            
            // Finds the Lagrange multiplier step 
            static void lagrangeMultiplierStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & grad=state.grad;
                X_Vector const & dx=state.dx;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                X_Vector & x=state.x;
                Y_Vector & dy=state.dy;

                // x_p_dx <- x + dx
                X_Vector x_p_dx(X::init(x));
                    X::copy(x,x_p_dx);
                    X::axpy(Real(1.),dx,x_p_dx);

                // grad_xpdx <- L(x+dx,y) = grad f(x+dx) + g'(x+dx)*y
                X_Vector grad_xpdx(X::init(x));
                    f.grad(x_p_dx,grad_xpdx);

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_xpdx_mult(X::init(grad));
                    f_mod.grad_mult(x_p_dx,grad_xpdx,grad_xpdx_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(dy));
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x+dx,y),0);
                XxY_Vector b0(XxY::init(x0));
                    X::copy(grad_xpdx_mult,b0.first);
                    X::scal(Real(-1.),b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // This is a somewhat unsatisfying hack.  Basically, many
                // of our operators like the preconditioner need to know
                // what the current iterate is.  Right now, we just have
                // a generic operator, which has no knowledge of this, so
                // we create references behind the scenes to link to the
                // current iterate.  Most of the time, this works fine,
                // except here were we're building our augmented system
                // at x+dx, but the preconditioner may and probably is linked
                // to x.  Hence, we're going to temporarily move where our
                // current iterate is for this solve and then move back.
                X_Vector x_save(X::init(x));
                    X::copy(x,x_save);
                X::copy(x_p_dx,x);

                // Solve the augmented system for the Lagrange multiplier step 
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    LagrangeMultiplierStepManipulator(state,fns),
                    x0 
                );

                // Restore our current iterate
                X::copy(x_save,x);

                // Copy out the Lagrange multiplier step
                Y::copy(x0.second,dy);
            }
            
            // Does a check on how far off the Lagrange multiplier is 
            static Real lagrangeMultiplierCheck(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & grad=state.grad; 
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                X_Vector & x=state.x;
                Y_Vector & dy=state.dy;

                // grad_x <- L(x,y) = grad f(x) + g'(x)*y
                X_Vector grad_x(X::init(x));
                    f.grad(x,grad_x);

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_x_mult(X::init(grad));
                    f_mod.grad_mult(x,grad_x,grad_x_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0(X::init(x),Y::init(dy));
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x,y),0);
                XxY_Vector b0(XxY::init(x0));
                    X::copy(grad_x_mult,b0.first);
                    X::scal(Real(-1.),b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the Lagrange multiplier step 
                Optizelle::gmres <Real,XXxYY> (
                    AugmentedSystem(state,fns,x),
                    b0,
                    Real(1.), // This will be overwritten by the manipulator
                    augsys_iter_max,
                    augsys_rst_freq,
                    PAugSys_l,
                    PAugSys_r,
                    LagrangeMultiplierStepManipulator(state,fns),
                    x0 
                );

                // Copy out the Lagrange multiplier step
                return sqrt(Y::innr(x0.second,x0.second));
            }

            // Computes the predicted reduction 
            static void predictedReduction(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                X_Vector const & W_gradpHdxn=state.W_gradpHdxn;
                X_Vector const & H_dxn=state.H_dxn;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & H_dxtuncorrected=state.H_dxtuncorrected;
                X_Vector const & dx_t_uncorrected=state.dx_t_uncorrected;
                Y_Vector const & gpxdxn_p_gx=state.gpxdxn_p_gx;
                Y_Vector const & dy=state.dy;
                Y_Vector const & g_x = state.g_x;
                Real const & rho=state.rho;
                Real const & norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                Real & pred=state.pred;
                
                // Find || g(x) ||
                Real norm_gx = sqrt(Y::innr(g_x,g_x));
                
                // Find the gradient modifications for step computation 
                X_Vector grad_step(X::init(grad));
                    f_mod.grad_step(x,grad,grad_step);
                
                // Add the Hessian modifications to H(x)dx_n
                X_Vector Hdxn_step(X::init(x));
                    f_mod.hessvec_step(x,dx_n,H_dxn,Hdxn_step);
                
                // Add the Hessian modifications to H(x)dx_t_uncorrected
                X_Vector H_dxtuncorrected_step(X::init(x));
                    f_mod.hessvec_step(x,dx_t_uncorrected,H_dxtuncorrected,
                        H_dxtuncorrected_step);

                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                pred = - X::innr(W_gradpHdxn,dx_t_uncorrected);

                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                //    - .5 < H dx_t_uncorrected,dx_t_uncorrected >
                pred-=Real(0.5)*X::innr(H_dxtuncorrected_step,dx_t_uncorrected);

                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                //    - .5 < H dx_t_uncorrected,dx_t_uncorrected >
                //    - < grad L(x,y), dx_n >
                pred -= X::innr(grad_step,dx_n);
                
                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                //    - .5 < H dx_t_uncorrected,dx_t_uncorrected >
                //    - < grad L(x,y), dx_n > - .5 < H dx_n,dx_n >
                pred -= Real(0.5) * X::innr(Hdxn_step,dx_n); 
                
                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                //    - .5 < H dx_t_uncorrected,dx_t_uncorrected >
                //    - < grad L(x,y), dx_n > - .5 < H dx_n,dx_n >
                //    - < dy , g'(x)dx_n+g(x) >
                pred -= Y::innr(dy,gpxdxn_p_gx);
                
                // pred <- - < W (grad L(x,y) + H dx_n), dx_t_uncorrected >
                //    - .5 < H dx_t_uncorrected,dx_t_uncorrected >
                //    - < grad L(x,y), dx_n > - .5 < H dx_n,dx_n >
                //    - < dy , g'(x)dx_n+g(x) >
                //    + rho ( || g(x) ||^2 - || g'(x)dx_n+g(x) ||^2 )
                pred += rho* (norm_gx*norm_gx - norm_gpxdxnpgx*norm_gpxdxnpgx);
            }
            
            // Computes the penalty parameter 
            static void penaltyParameter(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Y_Vector const & g_x=state.g_x;
                Real const & pred=state.pred;
                Real const & norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                Real const & rho_old=state.rho_old;
                Real const & rho_bar=state.rho_bar;
                Real & rho=state.rho;
               
                // norm_gx <- || g(x) ||
                Real const & norm_gx=sqrt(Y::innr(g_x,g_x));

                // If the predicted reduction is small, update the penalty
                // parameter.  Make sure we actually have a positive predicted
                // correction before we attempt this.
                if( pred > 0 && 
                    pred < (rho_old/Real(2.))
                        * (norm_gx*norm_gx - norm_gpxdxnpgx*norm_gpxdxnpgx) 
                ) {
                    rho = -Real(2.) * pred
                        / (norm_gx*norm_gx - norm_gpxdxnpgx*norm_gpxdxnpgx) 
                        + Real(2.) * rho_old + rho_bar;
                }
            }

            // Computes the residual predicted reduction 
            static void residualPredictedReduction(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Y_Vector const & dy=state.dy;
                Y_Vector const & gpxdxn_p_gx=state.gpxdxn_p_gx;
                Y_Vector const & gpxdxt=state.gpxdxt;
                Real const & rho=state.rho;
                Real & rpred=state.rpred;

                // rpred <- - < dy, g'(x)dx_t>
                rpred = -Y::innr(dy,gpxdxt);

                // rpred <- - < dy, g'(x)dx_t> - rho || g'(x)dx_t ||^2
                rpred -= rho * Y::innr(gpxdxt,gpxdxt);

                // rpred <- - < dy, g'(x)dx_t> - rho || g'(x)dx_t ||^2
                //     - 2 rho < g'(x)dx_t, g'(x) dx_n + g(x) >
                rpred -= Real(2.) * rho * Y::innr(gpxdxt,gpxdxn_p_gx);
            }
            
            // Checks whether we accept or reject a step
            static bool checkStep(
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create shortcuts to some elements in the state
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                Y_Vector const & dy=state.dy;
                Real const & eta1=state.eta1;
                Real const & eta2=state.eta2;
                Real const & f_x=state.f_x;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Y_Vector & y=state.y;
                Real & delta=state.delta;
                Real & ared=state.ared;
                Real & pred=state.pred;
                Real & f_xpdx=state.f_xpdx;
                
                // Allocate memory for temporaries that we need
                X_Vector x_p_dx(X::init(x));

                // Determine x+dx 
                X::copy(dx,x_p_dx);
                X::axpy(Real(1.),x,x_p_dx);

                // Save the old Lagrange multiplier
                Y_Vector y_old(Y::init(y));
                    Y::copy(y,y_old);

                // Determine y + dy
                Y::axpy(Real(1.),dy,y);

                // Determine the merit function at x and x+dx
                Real merit_x = f_mod.merit(x,f_x);
                f_xpdx = f(x_p_dx);
                Real merit_xpdx = f_mod.merit(x_p_dx,f_xpdx);

                // Restore the old Lagrange multiplier
                Y::copy(y_old,y);

                // norm_dx = || dx ||
                Real norm_dx = sqrt(X::innr(dx,dx));

                // Determine the actual reduction
                ared = merit_x - merit_xpdx;
                
                // Add a safety check in case we don't actually minimize the TR
                // subproblem correctly. This could happen for a variety of
                // reasons.  Most notably, if we do not correctly calculate the
                // Hessian approximation, we could have a nonsymmetric 
                // approximation.  In that case, truncated-CG will exit, but 
                // has an undefined result.  In the case that the actual 
                // reduction also increases, rho could have an extraneous 
                // positive value.  Hence, we require an extra check.
                if(pred < Real(0.)){
                    delta = norm_dx/Real(2.);
                    return false;
                }

                // Update the trust region radius and return whether or not we
                // accept the step
                if(ared >= eta2*pred){
                    // Increase the size of the trust-region if the Krylov
                    // solver reached the boundary.
                    if( krylov_stop==KrylovStop::NegativeCurvature ||
                        krylov_stop==KrylovStop::TrustRegionViolated
                    ) 
                        delta *= Real(2.);
                    return true;
                } else if(ared >= eta1*pred && ared < eta2*pred)
                    return true;
                else {
                    delta = norm_dx/Real(2.);
                    return false;
                }
            }

            // Finds the trust-region step
            static void getStep(
                Messaging const & msg,
                StateManipulator <EqualityConstrained <Real,XX,YY> > const &
                    smanip,
                typename Functions::t const & fns,
                typename State::t & state
            ){
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);
                X_Vector const & x=state.x;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & dx_t=state.dx_t;
                X_Vector const & dx_tcp_uncorrected
                    =state.dx_tcp_uncorrected;
                Real const & xi_4=state.xi_4;
                Real const & eta0=state.eta0;
                Real const & eps_dx=state.eps_dx;
                Real const & norm_dxtyp=state.norm_dxtyp;
                Real const & rho_old=state.rho_old;
                Natural const & history_reset=state.history_reset;
                X_Vector & dx=state.dx;
                X_Vector & dx_t_uncorrected=state.dx_t_uncorrected;
                X_Vector & H_dxn=state.H_dxn;
                X_Vector & H_dxtuncorrected=state.H_dxtuncorrected;
                Y_Vector & g_x=state.g_x;
                Y_Vector & gpxdxn_p_gx=state.gpxdxn_p_gx;
                Y_Vector & gpxdxt=state.gpxdxt;
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                Real & norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                Real & xi_qn=state.xi_qn;
                Real & xi_pg=state.xi_pg;
                Real & xi_proj=state.xi_proj;
                Real & xi_tang=state.xi_tang;
                Real & xi_lmh=state.xi_lmh;
                Real & pred=state.pred;
                Real & rpred=state.rpred;
                Real & rho=state.rho;
                Real & alpha0=state.alpha0;
                Natural & rejected_trustregion=state.rejected_trustregion;

                // Create a single temporary vector
                X_Vector x_tmp1(X::init(x));

                // Continue to look for a step until our actual vs. predicted
                // reduction is good.
                rejected_trustregion=0;
                while(true) {
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeGetStep);

                    // Continue to look for a step until the inaccuracy
                    // in the normal and tangential steps are acceptable.
                    // The iterate i alternates between trying the Cauchy
                    // point in the tangential step and recomputing the
                    // tangential step.  If we're on an even iteration, we
                    // have a freshly computed tangential step.  If we're on
                    // an odd iteration, we're attempting to use the Cauchy
                    // point.
                    for(int i=0;;i++) { 

                        // Compute a brand new step
                        if(i%2 == 0) {
                            // Find the quasi-Normal step
                            quasinormalStep(fns,state);

                            // Find g'(x) dx_n + g(x)
                            g.p(x,dx_n,gpxdxn_p_gx);
                            Y::axpy(Real(1.),g_x,gpxdxn_p_gx);

                            // Find || g'(x) dx_n + g(x) ||
                            norm_gpxdxnpgx
                                = sqrt(Y::innr(gpxdxn_p_gx,gpxdxn_p_gx));

                            // Find H dx_n
                            f.hessvec(x,dx_n,H_dxn);
                            
                            // Find W (g + H dxn) 
                            projectedGradLagrangianPlusHdxn(fns,state);

                            // Find the uncorrected tangential step
                            tangentialSubProblem(fns,state);
                        }
                    
                        // Find H dx_t_uncorrected
                        f.hessvec(x,dx_t_uncorrected,H_dxtuncorrected);

                        // Continue to correct the tangential step until one
                        // comes back as valid.  Sometimes, we can get a
                        // negative predicted reduction due to odd numerical
                        // errors or a bad Hessian.  In this case, no amount
                        // of refining can fix the problem, so just exit and
                        // let the checkStep code adjust the trust-region.
                        rpred=Real(1.);
                        pred=Real(0.);
                        Real xi_tang0 = xi_tang;
                        for( ;
                            pred>=Real(0.) && fabs(rpred)>eta0*pred
                                && xi_tang>std::numeric_limits<Real>::epsilon();
                            xi_tang=xi_tang*Real(.001)
                        ) {

                            // Correct the tangential step
                            tangentialStep(fns,state);

                            // Find the primal step
                            X::copy(dx_n,dx);
                            X::axpy(Real(1.),dx_t,dx);

                            // Find g'(x)dx_t
                            g.p(x,dx_t,gpxdxt);

                            // Find the Lagrange multiplier step
                            lagrangeMultiplierStep(fns,state);

                            // Find the predicted reduction
                            rho = rho_old;
                            predictedReduction(fns,state);

                            // Update the penalty paramter
                            penaltyParameter(fns,state);

                            // Find the predicted reduction based on the new
                            // penalty parameter
                            predictedReduction(fns,state);

                            // Find the residual predicted reduction
                            residualPredictedReduction(fns,state);
                        }
                        xi_tang = xi_tang0;

                        // Check if ||dx_t_uncorrected||<=xi_4 || dx_n + dx_t||.
                        // In this case, the inexactness is acceptable and we
                        // can exit.
                        if( X::innr(dx_t_uncorrected,dx_t_uncorrected) <=
                                xi_4*xi_4 * X::innr(dx,dx)
                        ) break;

                        // If the inexactness isn't acceptable, try the Cauchy
                        // point before recomputing everything 
                        if(i % 2==0)
                            X::copy(dx_tcp_uncorrected,dx_t_uncorrected);

                        // If the Cauchy point didn't work, then tighten the
                        // tolerances and try again
                        else {
                            xi_qn /= Real(10.);
                            xi_pg /= Real(10.);
                            xi_proj /= Real(10.);
                            xi_tang /= Real(10.);
                            xi_lmh /= Real(10.);

                            // If any tolerance hits essentially zero, set
                            // the step to zero and allow the algorithm to
                            // exit
                            if(xi_qn < std::numeric_limits <Real>::epsilon() ||
                                xi_pg < std::numeric_limits <Real>::epsilon() ||
                                xi_proj<std::numeric_limits <Real>::epsilon() ||
                                xi_tang<std::numeric_limits <Real>::epsilon() ||
                                xi_lmh <std::numeric_limits <Real>::epsilon()
                            ){
                                X::zero(dx);
                                pred=Real(0.);
                            }
                        }
                    }

                    // In an interior point method, we may truncate the step
                    // and this information is communicated through the
                    // base linesearch parameter, alpha0.
                    alpha0 = Real(1.);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::BeforeActualVersusPredicted);

                    // If need be, shorten the step
                    X::scal(alpha0,dx);

                    // If we shorten our step, update our Lagrange multiplier
                    // step
                    if(alpha0 < Real(1.))
                        lagrangeMultiplierStep(fns,state);
                    
                    // Check whether the step is good
                    if(checkStep(fns,state))
                        break;
                    else
                        rejected_trustregion++;
                    
                    // If the number of rejected steps is above the
                    // history_reset threshold, destroy the quasi-Newton
                    // information
                    if(rejected_trustregion > history_reset){
                        oldY.clear();
                        oldS.clear();
                    }

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::AfterRejectedTrustRegion);

                    // Alternatively, check if the step becomes so small
                    // that we're not making progress.  In this case, break
                    // and allow the stopping conditions to terminate
                    // optimization.  We use a zero length step so that we
                    // do not modify the current iterate.
                    Real norm_dx = sqrt(X::innr(dx,dx));
                    if(norm_dx < eps_dx*norm_dxtyp) {
                        X::scal(Real(0.),dx);
                        break;
                    }
                } 
            }
            
            // Adjust the stopping conditions unless
            // || g(x) || <  eps_constr || g(x_0) ||
            static void adjustStoppingConditions(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Y_Vector const & g_x = state.g_x;
                Real const & eps_constr=state.eps_constr;
                Real const & norm_gxtyp=state.norm_gxtyp; 
                StoppingCondition::t & opt_stop=state.opt_stop;
                
                // Prevent convergence unless the infeasibility is small. 
                Real norm_gx=sqrt(Y::innr(g_x,g_x));
                if( opt_stop==StoppingCondition::RelativeGradientSmall &&
                    !(norm_gx < eps_constr*norm_gxtyp) 
                )
                    opt_stop=StoppingCondition::NotConverged;
            }

            
            // This adds the composite-step method through use of a state
            // manipulator.
            template <typename ProblemClass>
            struct CompositeStepManipulator
                : public StateManipulator <ProblemClass>
            {
            private:
                // A reference to the user-defined state manipulator
                StateManipulator<ProblemClass> const & smanip;

                // A reference to the messaging object
                Messaging const & msg;

            public:
                CompositeStepManipulator(
                    StateManipulator <ProblemClass> const & smanip_,
                    Messaging const & msg_
                ) : smanip(smanip_), msg(msg_) {}


                // Application
                void operator () (
                    typename ProblemClass::Functions::t const & fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t const & loc
                ) const {
                    // Call the user define manipulator
                    smanip(fns_,state_,loc);

                    // Dynamically cast the incoming state and fns to the
                    // to work with the equality constrained spaces.  In theory,
                    // this should always work since we're doing this trickery
                    // internally.  Basically, this is required since we're
                    // inserting into the unconstrained constrained code.
                    // Within this code, the state manipulator is hard coded to 
                    // use the state for the unconstrained problem even though
                    // this state is really an equality constrained state when
                    // called using the routines below.
                    typename Functions::t const & fns
                        =dynamic_cast <typename Functions::t const &> (fns_);
                    typename State::t & state 
                        =dynamic_cast <typename State::t &> (state_);
                
                    // Create some shortcuts
                    ScalarValuedFunctionModifications <Real,XX> const & f_mod
                        = *(fns.f_mod);
                    VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);
                    X_Vector const & x=state.x;
                    Y_Vector const & dy=state.dy;
                    Real const & rho = state.rho;
                    Natural const & krylov_iter_max = state.krylov_iter_max;
                    AlgorithmClass::t& algorithm_class=state.algorithm_class;
                    X_Vector & grad=state.grad;
                    Y_Vector & y=state.y;
                    Y_Vector & g_x=state.g_x;
                    Real & norm_gxtyp = state.norm_gxtyp;
                    Real & rho_old = state.rho_old;
                    Real & norm_gradtyp = state.norm_gradtyp;
                    Natural & krylov_orthog_max = state.krylov_orthog_max;

                    switch(loc){
                    case OptimizationLocation::BeforeInitialFuncAndGrad:
                        // Make sure the algorithm uses the composite step
                        // routines to find the new step.
                        algorithm_class=AlgorithmClass::UserDefined;

                        // Make sure that we do full orthogonalization in our
                        // truncated Krylov method
                        krylov_orthog_max = krylov_iter_max;
                        break;
                
                    case OptimizationLocation::AfterInitialFuncAndGrad: {
                        // Make sure we properly cache g(x) and its norm
                        // on initialization.  
                        g(x,g_x);
                        norm_gxtyp = sqrt(Y::innr(g_x,g_x));

                        // Find the initial Lagrange multiplier and then update
                        // the gradient and merit function.
                        lagrangeMultiplier(fns,state);

                        // In addition, update the norm of gradient and
                        // typical gradient since we've modified the Lagrange
                        // multiplier
                        X_Vector grad_stop(X::init(grad));
                            f_mod.grad_stop(x,grad,grad_stop);
                        norm_gradtyp=sqrt(X::innr(grad_stop,grad_stop));

                        // Prime the previous penalty parameter
                        rho_old = rho;
                        break;

                    } case OptimizationLocation::GetStep: {
                        // Get a manipulator compatible with equality
                        // constrained code
                        ConversionManipulator
                            <ProblemClass,EqualityConstrained<Real,XX,YY> >
                            cmanip(smanip);

                        // Find the steps in both the primal and dual directions
                        getStep(msg,cmanip,fns,state);
                        break;

                    } case OptimizationLocation::BeforeStep:
                        // Save the new penalty parameter
                        rho_old = rho; 

                        // Make sure to take the step in the dual variable
                        Y::axpy(Real(1.),dy,y);
                        break;

                    case OptimizationLocation::AfterGradient: {
                        // In an interior point method, we may have modified
                        // our interior point parameter, which changes the
                        // gradient.  This necessitates a new Lagrange 
                        // multiplier computation. 
                        lagrangeMultiplier(fns,state);
                        break;

                    } case OptimizationLocation::AfterStepBeforeGradient:
                        // Make sure we update our cached value of g(x) 
                        g(x,g_x);
                        break;

                    case OptimizationLocation::EndOfOptimizationIteration:
                        // Make sure we don't exit until the norm of the
                        // constraints is small as well.
                        adjustStoppingConditions(fns,state);
                        break;

                    default:
                        break;
                    }
                }
            };

            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                StateManipulator <EqualityConstrained <Real,XX,YY> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                StateManipulator <EqualityConstrained <Real,XX,YY> > const &
                    smanip,
                typename Functions::t & fns,
                typename State::t & state
            ){
                
                // Adds the output pieces to the state manipulator 
                DiagnosticManipulator <EqualityConstrained <Real,XX,YY> >
                    dmanip(smanip,msg);

                // Add the composite step pieces to the state manipulator
                CompositeStepManipulator <EqualityConstrained <Real,XX,YY> >
                    csmanip(dmanip,msg);

                // Insures that we can interact with unconstrained code
                ConversionManipulator
                    <EqualityConstrained<Real,XX,YY>,Unconstrained <Real,XX> >
                    cmanip(csmanip);
                
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);
                
                // Check the inputs to the optimization
                State::check(msg,state);

                // Minimize the problem
                Unconstrained <Real,XX>::Algorithms
                    ::getMin_(msg,cmanip,fns,state);
            }
        };
    };
        
    // Routines that manipulate and support problems of the form
    // 
    // min_{x \in X} f(x) st h(x) >=_K 0
    //
    // where f : X -> R and h : X -> Z
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class ZZ
    > 
    struct InequalityConstrained {
        // Disallow constructors
        NO_CONSTRUCTORS(InequalityConstrained);

        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State);

            // Internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t);

                // Inequality multiplier (dual variable or Lagrange multiplier)
                Z_Vector z;
                
                // Step in the inequality multiplier 
                Z_Vector dz;

                // The inequality constraint evaluated at x.  In theory,
                // we can always just evaluate this when we need it.  However,
                // we do require its computation both in the gradient as well
                // as Hessian calculations.  More specifically, when computing
                // with SDP constraints, we require a factorization of this
                // quantity.  By caching it, we have the ability to cache the
                // factorization.
                Z_Vector h_x;

                // Interior point parameter
                Real mu;

                // Current interior point estimate
                Real mu_est;

                // Typical value for mu.  Generally, the first estimated
                // value for mu.
                Real mu_typ;

                // Relative stopping criteria for the interior point parameter
                Real eps_mu;

                // The amount that we reduce the interior point parameter by
                // everytime we approach the central path
                Real sigma;

                // How close we move to the boundary during a single step
                Real gamma;

                // Type of interior point method
                InteriorPointMethod::t ipm;

                // Centrality strategy
                CentralityStrategy::t cstrat;

                // Initialization constructors
                t(X_Vector const & x_,Z_Vector const & z_) :
                    Unconstrained <Real,XX>::State::t(x_),
                    z(Z::init(z_)),
                    dz(Z::init(z_)),
                    h_x(Z::init(z_)),
                    mu(std::numeric_limits<Real>::quiet_NaN()),
                    mu_est(std::numeric_limits<Real>::quiet_NaN()),
                    mu_typ(std::numeric_limits<Real>::quiet_NaN()),
                    eps_mu(1e-8),
                    sigma(0.5),
                    gamma(0.95),
                    ipm(InteriorPointMethod::PrimalDual),
                    cstrat(CentralityStrategy::Constant)
                {
                    Z::copy(z_,z);
                }
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };
                
            // Check that we have a valid set of parameters.  
            static void check_(Messaging const & msg,t const & state) {
                // Use this to build an error message
                std::stringstream ss;
                
                // Check that the interior point parameter is positive 
                if(state.mu <= Real(0.)) 
                    ss << "The interior point parameter must be positive: " 
                        "mu = " << state.mu;
                
                // Check that the interior point parameter estimate is positive 
                if(state.mu_est <= Real(0.)) 
                    ss << "The interior point parameter estimate must be "
                        "positive: mu_est = " << state.mu_est;

                // Check that the typical interior point parameter is positive 
                if(state.mu_typ <= Real(0.)) 
                    ss << "The typical interior point parameter must be "
                        "positive:  mu_typ = " << state.mu_typ;

                // Check that the interior point stopping tolerance is positive 
                else if(state.eps_mu <= Real(0.)) 
                    ss << "The interior point stopping tolerance must be "
                        "positive: eps_mu = " << state.eps_mu;

                // Check that the reduction in the interior point parameter
                // is between 0 and 1.
                else if(state.sigma <= Real(0.) || state.sigma >= Real(1.)) 
                    ss << "The reduction in the interior point parameter "
                        "must be between 0 and 1: sigma = " << state.sigma;

                // Check that the fraction to the boundary is between 0 and 1. 
                else if(state.gamma <= Real(0.) || state.gamma >= Real(1.)) 
                    ss << "The fraction to the boundary must be between " 
                        "0 and 1: gamma= " << state.gamma;

                // If there's an error, print it
                if(ss.str()!="") msg.error(ss.str());
            }
            static void check(Messaging const & msg,t const & state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                InequalityConstrained <Real,XX,ZZ>::State::check_(msg,state);
            }
        };
        // Utilities for restarting the optimization
        struct Restart {
            // Disallow constructors
            NO_CONSTRUCTORS(Restart);
       
            // Holds restart information
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <Natural> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Z_Vector> > Z_Vectors;
            
            // Checks whether we have a valid real label.
            struct is_real :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_real()(name) ||
                        name == "mu" ||
                        name == "mu_est" ||
                        name == "mu_typ" ||
                        name == "eps_mu" ||
                        name == "sigma" ||
                        name == "gamma" 
                    )
                        return true;
                    else
                        return false;
                    }
            };
            
            // Checks whether we have a valid natural number label.
            struct is_nat :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                        ::is_nat()(name)
                    )
                        return true;
                    else
                        return false;
                }
            };
           
            // Checks whether we have a valid parameter label.
            struct is_param :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_param()(name) ||
                        name == "ipm" ||
                        name == "cstrat"
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid variable label
            struct is_x : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_x()(name) 
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid inequality multiplier label
            struct is_z : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( name == "z" ||
                        name == "dz" ||
                        name == "h_x"
                    )
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have valid labels
            static void checkLabels(
                Messaging const & msg,
                Reals const & reals,
                Nats const & nats,
                Params const & params,
                X_Vectors const & xs,
                Z_Vectors const & zs
            ) {
                Utility::checkLabels(msg,is_real(),reals.first," real name: ");
                Utility::checkLabels(msg,is_nat(),nats.first," natural name: ");
                Utility::checkLabels(
                    msg,is_param(),params.first," paramater name: ");
                Utility::checkLabels(msg,is_x(),xs.first," variable name: ");
                Utility::checkLabels(
                    msg,is_z(),zs.first," inequality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string const &,std::string const &,std::string>
            {
                std::string operator() (
                    std::string const & label,
                    std::string const & val
                ) const {

                    // Create a base message
                    const std::string base
                        ="During serialization, found an invalid ";

                    // Used to build the message 
                    std::stringstream ss;

                    // Check the unconstrained parameters
                    if(typename Unconstrained <Real,XX>
                        ::Restart::is_param()(label)
                    ) {
                        ss << typename Unconstrained <Real,XX>::Restart
                            ::checkParamVal()(label,val);

                    // Check the interior point method type
                    } else if(label=="ipm") {
                        if(!InteriorPointMethod::is_valid()(val))
                            ss << base << "interior point method type: " << val;
                    } else if(label=="cstrat") {
                        if(!CentralityStrategy::is_valid()(val))
                            ss << base << "centrality strategy: " << val;
                    }
                    return ss.str();
                }
            };
            
            // Copy out the inequality multipliers 
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs,
                Z_Vectors & zs
            ) {
                zs.first.emplace_back("z");
                zs.second.emplace_back(std::move(state.z));
                zs.first.emplace_back("dz");
                zs.second.emplace_back(std::move(state.dz));
                zs.first.emplace_back("h_x");
                zs.second.emplace_back(std::move(state.h_x));
            }
            
            // Copy out the scalar information
            static void stateToScalars(
                typename State::t & state,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy in all the real numbers
                reals.first.emplace_back("mu");
                reals.second.emplace_back(state.mu);
                reals.first.emplace_back("mu_est");
                reals.second.emplace_back(state.mu_est);
                reals.first.emplace_back("mu_typ");
                reals.second.emplace_back(state.mu_typ);
                reals.first.emplace_back("eps_mu");
                reals.second.emplace_back(state.eps_mu);
                reals.first.emplace_back("sigma");
                reals.second.emplace_back(state.sigma);
                reals.first.emplace_back("gamma");
                reals.second.emplace_back(state.gamma);

                // Copy in all of the parameters
                params.first.emplace_back("ipm");
                params.second.emplace_back(
                    InteriorPointMethod::to_string(state.ipm));
                params.first.emplace_back("cstrat");
                params.second.emplace_back(
                    CentralityStrategy::to_string(state.cstrat));
            }
            
            // Copy in inequality multipliers 
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs,
                Z_Vectors & zs
            ) {
                typename std::list <Z_Vector>::iterator z
                    =zs.second.begin();
                for(typename std::list <std::string>::iterator name
                        =zs.first.begin();
                    name!=zs.first.end();
                ) {
                    // Make a copy of the current iterators.  We use these
                    // to remove elements
                    typename std::list <std::string>::iterator name0 = name;
                    typename std::list <Z_Vector>::iterator z0 = z;

                    // Increment our primary iterators 
                    name++; z++;

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name0=="z")
                        state.z = std::move(*z0);
                    else if(*name0=="dz")
                        state.dz = std::move(*z0);
                    else if(*name0=="h_x")
                        state.h_x = std::move(*z0);
                }
            }
            
            // Copy in the scalar information
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Nats & nats,
                Params & params
            ) { 
                // Copy in any reals 
                typename std::list <Real>::iterator real=reals.second.begin();
                for(std::list <std::string>::iterator name=reals.first.begin();
                    name!=reals.first.end();
                    name++,real++
                ){
                    if(*name=="mu") state.mu=*real;
                    else if(*name=="mu_est") state.mu_est=*real;
                    else if(*name=="mu_typ") state.mu_typ=*real;
                    else if(*name=="eps_mu") state.eps_mu=*real;
                    else if(*name=="sigma") state.sigma=*real;
                    else if(*name=="gamma") state.gamma=*real;
                } 
                    
                // Next, copy in any parameters 
                std::list <std::string>::iterator param=params.second.begin();
                for(std::list <std::string>::iterator name=params.first.begin();
                    name!=params.first.end();
                    name++,param++
                ){
                    if(*name=="ipm")
                        state.ipm=InteriorPointMethod::from_string(*param);
                    else if(*name=="cstrat")
                        state.cstrat=CentralityStrategy::from_string(*param);
                }
            }

            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Z_Vectors & zs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>
                    ::Restart::stateToVectors(state,xs);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::stateToVectors(state,xs,zs);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::stateToScalars(state,reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::stateToScalars(state,reals,nats,params);
            }
            
            // Capture data from structures controlled by the user.  
            static void capture(
                Messaging const & msg,
                typename State::t & state,
                X_Vectors & xs,
                Z_Vectors & zs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,zs);

                // Check the strings used to represent parameters
                Utility::checkParams(msg,checkParamVal(),params);

                // Copy in the variables 
                Unconstrained <Real,XX>
                    ::Restart::vectorsToState(state,xs);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::vectorsToState(state,xs,zs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::scalarsToState(state,reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::scalarsToState(state,reals,nats,params);

                // Check that we have a valid state 
                State::check(msg,state);
            }
        };
        
        // All the functions required by an optimization algorithm.  Note, this
        // routine owns the memory for these operations.  
        struct Functions {
            // Disallow constructors
            NO_CONSTRUCTORS(Functions);

            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t);

                // Inequality constraints 
                std::unique_ptr <VectorValuedFunction <Real,XX,ZZ> > h;
                
                // Initialize all of the pointers to null
                t() : Unconstrained <Real,XX>::Functions::t(), h(nullptr) {}
            };

            struct InequalityModifications
                : public Optizelle::ScalarValuedFunctionModifications <Real,XX>
            {
            public:
                // Disallow constructors
                NO_COPY_ASSIGNMENT(InequalityModifications); 

            private:
                // Underlying modification.  This takes control of the memory
                std::unique_ptr <
                    Optizelle::ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Inequality constraint.
                Optizelle::VectorValuedFunction <Real,XX,ZZ> const & h;
                
                // Inequality Lagrange multiplier
                Z_Vector const & z;

                // Interior point parameter
                Real const & mu;

                // Inequality constraint evaluated at x
                Z_Vector const & h_x;
                
                // Some workspace for the below functions
                mutable X_Vector grad_tmp;
                mutable X_Vector hess_mod; 
                mutable X_Vector x_tmp1;
                mutable Z_Vector z_tmp1;
                mutable Z_Vector z_tmp2;
                
                // Variables used for caching.  The boolean values denote
                // whether or not we've started caching yet.
                mutable std::pair <bool,X_Vector> x_merit;
                mutable Z_Vector hx_merit;
                mutable std::pair <bool,X_Vector> x_lag;
                mutable std::pair <bool,Z_Vector> z_lag;
                mutable std::pair <bool,X_Vector> x_schur;
                mutable std::pair <bool,Z_Vector> z_schur;
                mutable X_Vector hpxsz;
                mutable X_Vector hpxs_invLhx_e;

                // Adds the Lagrangian pieces to the gradient
                void grad_lag(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_lag
                ) const {
                    // grad_lag <- grad f(x)
                    X::copy(grad,grad_lag);
                    
                    // If relative error between the current and cached values
                    // is large, compute anew.
                    if( rel_err_cached <Real,XX> (x,x_lag)
                            >= std::numeric_limits <Real>::epsilon()*1e1 ||
                        rel_err_cached <Real,ZZ> (z,z_lag)
                            >= std::numeric_limits <Real>::epsilon()*1e1 
                    ) {
                        // hpxsz <- h'(x)* z 
                        h.ps(x,z,hpxsz);

                        // Cache the values
                        x_lag.first=true;
                        X::copy(x,x_lag.second);
                        z_lag.first=true;
                        Z::copy(z,z_lag.second);
                    }

                    // grad_lag <- grad f(x) - h'(x)*z
                    X::axpy(-Real(1.0),hpxsz,grad_lag);
                }
                
                // Adds the Schur complement pieces to the gradient
                void grad_schur(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_schur
                ) const {
                    // grad_schur <- grad f(x)
                    X::copy(grad,grad_schur);
                    
                    // If relative error between the current and cached values
                    // is large, compute anew.
                    if( rel_err_cached <Real,XX> (x,x_schur)
                            >= std::numeric_limits <Real>::epsilon()*1e1 ||
                        rel_err_cached <Real,ZZ> (z,z_schur)
                            >= std::numeric_limits <Real>::epsilon()*1e1
                    ) {
                        // z_tmp1 <- e
                        Z::id(z_tmp1);

                        // z_tmp2 <- inv(L(h(x))) e 
                        Z::linv(h_x,z_tmp1,z_tmp2);

                        // hpxs_invLhx_e <- h'(x)* (inv(L(h(x))) e)
                        h.ps(x,z_tmp2,hpxs_invLhx_e);
                        
                        // Cache the values
                        x_schur.first=true;
                        X::copy(x,x_schur.second);
                        z_schur.first=true;
                        Z::copy(z,z_schur.second);
                    }

                    // grad_schur<- grad f(x) - mu h'(x)* (inv(L(h(x))) e)
                    X::axpy(-mu,hpxs_invLhx_e,grad_schur);
                }
            public:
                InequalityModifications(
                    typename State::t const & state,
                    typename Functions::t & fns
                ) : f_mod(std::move(fns.f_mod)),
                    h(*(fns.h)),
                    z(state.z),
                    mu(state.mu),
                    h_x(state.h_x),
                    grad_tmp(X::init(state.x)),
                    hess_mod(X::init(state.x)),
                    x_tmp1(X::init(state.x)),
                    z_tmp1(Z::init(state.z)),
                    z_tmp2(Z::init(state.z)),
                    x_merit(false,X::init(state.x)),
                    hx_merit(Z::init(state.z)),
                    x_lag(false,X::init(state.x)),
                    z_lag(false,Z::init(state.z)),
                    x_schur(false,X::init(state.x)),
                    z_schur(false,Z::init(state.z)),
                    hpxsz(X::init(state.x)),
                    hpxs_invLhx_e(X::init(state.x))
                {}

                // Merit function additions to the objective
                virtual Real merit(X_Vector const & x,Real const & f_x) const {
                    // Do the underlying modification of the objective
                    Real merit_x = f_mod->merit(x,f_x);
                    
                    // If we've not started caching or the relative error
                    // is large, compute anew.
                    if( rel_err_cached <Real,XX> (x,x_merit)
                            >= std::numeric_limits <Real>::epsilon()*1e1
                    ) {
                        // hx_merit <- h(x)
                        h(x,hx_merit);
                        
                        // Cache the values
                        x_merit.first=true;
                        X::copy(x,x_merit.second);
                    }

                    // Return merit(x) - mu barr(h(x))
                    return merit_x - mu * Z::barr(hx_merit); 
                }

                // Stopping condition modification of the gradient
                virtual void grad_stop(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_stop
                ) const {
                    f_mod->grad_stop(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_stop);
                }

                // Diagnostic modification of the gradient
                virtual void grad_diag(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_diag
                ) const {
                    f_mod->grad_diag(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_diag);
                }

                // Modification of the gradient when finding a trial step
                virtual void grad_step(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_step
                ) const {
                    f_mod->grad_step(x,grad,grad_tmp);
                    grad_schur(x,grad_tmp,grad_step);
                }

                // Modification of the gradient for a quasi-Newton method 
                virtual void grad_quasi(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_quasi
                ) const {
                    f_mod->grad_quasi(x,grad,grad_quasi);
                }

                // Modification of the gradient when solving for the equality
                // multiplier
                virtual void grad_mult(
                    X_Vector const & x,
                    X_Vector const & grad,
                    X_Vector & grad_mult
                ) const {
                    f_mod->grad_mult(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_mult);
                }

                // Modification of the Hessian-vector product when finding a
                // trial step
                virtual void hessvec_step(
                    X_Vector const & x,
                    X_Vector const & dx,
                    X_Vector const & H_dx,
                    X_Vector & Hdx_step 
                ) const {

                    // Modify the Hessian-vector product
                    f_mod->hessvec_step(x,dx,H_dx,Hdx_step);

                    // z_tmp1 <- h'(x) dx
                    h.p(x,dx,z_tmp1);

                    // z_tmp2 <- h'(x) dx o z
                    Z::prod(z_tmp1,z,z_tmp2);

                    // linv_hx_hpx_prod_z <- inv(L(h(x))) (h'(x) dx o z) 
                    Z::linv(h_x,z_tmp2,z_tmp1);

                    // hess_mod <- h'(x)* (inv(L(h(x))) (h'(x) dx o z))
                    h.ps(x,z_tmp1,hess_mod);

                    // H_dx 
                    //  = hess f(x) dx + h'(x)* (inv(L(h(x))) (h'(x) dx o z))
                    X::axpy(Real(1.),hess_mod,Hdx_step);
                }
            };

            // Check that all the functions are defined
            static void check(Messaging const & msg,t const & fns) {

                // Check the unconstrained pieces
                Unconstrained <Real,XX>::Functions::check(msg,fns);
                
                // Check that the inequality constraints exist 
                if(fns.h.get()==nullptr)
                    msg.error("Missing the inequality constraint definition.");
            }

            // Initialize any missing functions for just inequality constrained 
            // optimization.
            static void init_(
                Messaging const & msg,
                typename State::t & state,
                t& fns
            ) {
                // Check that all functions are defined 
                check(msg,fns);

                // Modify the objective 
                fns.f_mod.reset(new InequalityModifications(state,fns));

                // Set the trust-region scaling
#if 0
                fns.TRS.reset(
                    new typename Algorithms::TrustRegionScaling(fns,state));
#else
                if(fns.TRS.get()==nullptr)
                    fns.TRS.reset(new typename Unconstrained <Real,XX>
                        ::Functions::Identity());
#endif
            }

            // Initialize any missing functions 
            static void init(
                Messaging const & msg,
                typename State::t & state,
                t& fns
            ) {
                Unconstrained <Real,XX>
                    ::Functions::init_(msg,state,fns);
                InequalityConstrained <Real,XX,ZZ>
                    ::Functions::init_(msg,state,fns);
            }
        };
        
        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
            // Disallow constructors
            NO_CONSTRUCTORS(Printer);

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                // Print out the current interior point parameter and
                // the estimate of the interior point parameter.
                out.emplace_back(Utility::atos("mu"));
                out.emplace_back(Utility::atos("mu_est"));
            }

            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer::getStateHeader_
                    (state,out);
            }

            // Gets the state information for output
            static void getState_(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                Real const & mu=state.mu; 
                Real const & mu_est=state.mu_est; 

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Interior point information
                out.emplace_back(Utility::atos(mu));
                out.emplace_back(Utility::atos(mu_est));

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=Utility::blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                bool const & noiter,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getState_(fns,state,blank,out);
            }
            
            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getKrylovHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) { }

            // Combines all of the Krylov information
            static void getKrylov(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getKrylov_(state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms);

            // An operator to reshape the trust-region radius and, hopefully,
            // keep us away from the boundary.
            struct TrustRegionScaling : public Operator <Real,XX,XX> {
            private:
                // The function h 
                VectorValuedFunction <Real,XX,ZZ> const & h;

                // The value h(x) 
                Z_Vector const & h_x;

                // The current iterate
                X_Vector const & x;

                // Work vectors
                mutable Z_Vector z_tmp1;
                mutable Z_Vector z_tmp2;

            public:
                TrustRegionScaling(
                    typename Functions::t const & fns,
                    typename State::t const & state
                ) : h(*(fns.h)),
                    h_x(state.h_x),
                    x(state.x),
                    z_tmp1(Z::init(state.z)),
                    z_tmp2(Z::init(state.z))
                {}

                void operator () (X_Vector const & dx,X_Vector & result) const{
                    // z_tmp1 <- h'(x) dx
                    h.p(x,dx,z_tmp1); 

                    // z_tmp2 <- inv L(h(x)) h'(x) dx
                    Z::linv(h_x,z_tmp1,z_tmp2);

                    // result <- h'(x)* inv L(h(x)) h'(x) dx
                    h.ps(x,z_tmp2,result);
                }
            };

            // Finds the new inequality Lagrange multiplier
            // z = inv L(h(x)) (-h'(x)dx o z + mu e)
            static void findInequalityMultiplierLinked(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Z_Vector const & h_x=state.h_x;
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                Real const & mu=state.mu;
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                Z_Vector & z=state.z;

                // z_tmp1 <- h'(x)dx
                Z_Vector z_tmp1(Z::init(z));
                    h.p(x,dx,z_tmp1);

                // z_tmp2 <- h'(x)dx o z
                Z_Vector z_tmp2(Z::init(z));
                    Z::prod(z_tmp1,z,z_tmp2);

                // z_tmp2 <- -h'(x)dx o z
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 <- e
                Z::id(z_tmp1);

                // z_tmp2 <- -h'(x)dx o z + mu e
                Z::axpy(mu,z_tmp1,z_tmp2);

                // z <- inv L(h(x)) (-h'(x)dx o z + mu e)
                Z::linv(h_x,z_tmp2,z);

                // Symmetrize the iterate
                Z::symm(z);
            }

            // Finds the new inequality Lagrange multiplier
            // z = mu inv L(h(x)) e 
            static void findInequalityMultiplierLogBarrier(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Z_Vector const & h_x=state.h_x;
                Real const & mu=state.mu;
                Z_Vector & z=state.z;

                // z_tmp1 <- e 
                Z_Vector z_tmp1(Z::init(z));
                    Z::id(z_tmp1);

                // z <- inv(L(h(x))) e 
                Z::linv(h_x,z_tmp1,z);

                // z <- mu inv(L(h(x))) e 
                Z::scal(mu,z);
                
                // Symmetrize the iterate 
                Z::symm(z);
            }

            // Finds the new inequality Lagrange multiplier step
            // dz = -z + inv L(h(x)) (-h'(x)dx o z + mu e)
            static void findInequalityMultiplierStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Z_Vector const & z=state.z;
                Z_Vector const & h_x=state.h_x;
                X_Vector const & x=state.x;
                X_Vector const & dx=state.dx;
                Real const & mu=state.mu;
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                Z_Vector & dz=state.dz;

                // z_tmp1 <- h'(x)dx
                Z_Vector z_tmp1(Z::init(z));
                h.p(x,dx,z_tmp1);

                // z_tmp2 <- h'(x)dx o z
                Z_Vector z_tmp2(Z::init(z));
                Z::prod(z_tmp1,z,z_tmp2);

                // z_tmp2 <- -h'(x)dx o z
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 <- e
                Z::id(z_tmp1);

                // z_tmp2 <- -h'(x)dx o z + mu e
                Z::axpy(mu,z_tmp1,z_tmp2);

                // dz <- inv L(h(x)) (-h'(x)dx o z + mu e)
                Z::linv(h_x,z_tmp2,dz);

                // dz <- -z + inv L(h(x)) (-h'(x)dx o z + mu e)
                Z::axpy(Real(-1.),z,dz);
                
                // Symmetrize the direction
                Z::symm(dz);
            }

            // Estimates the interior point parameter with the formula
            // mu = <z,h(x)>/m
            static void estimateInteriorPointParameter(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Z_Vector const & z=state.z;
                Z_Vector const & h_x=state.h_x;
                Real & mu_est=state.mu_est;

                // Determine the scaling factor for the interior-
                // point parameter estimate
                Z_Vector z_tmp(Z::init(z));
                Z::id(z_tmp);
                Real m = Z::innr(z_tmp,z_tmp);

                // Estimate the interior-point parameter
                mu_est = Z::innr(z,h_x) / m;
            }

            // Find interior point parameter
            static void findInteriorPointParameter(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                Real const & mu_est=state.mu_est;
                Real const & mu_typ=state.mu_typ;
                Real const & sigma=state.sigma;
                Real const & eps_mu=state.eps_mu;
                CentralityStrategy::t const & cstrat=state.cstrat;
                Real const & norm_gradtyp=state.norm_gradtyp;
                Real const & eps_grad=state.eps_grad;
                Natural const & iter=state.iter;
                Real const & f_x=state.f_x;
                Real & mu=state.mu;
               
                // If we satisfy the stopping criteria, stop trying to
                // reduce the interior point parameter
                if(mu_est <= mu_typ*eps_mu) {
                    mu=mu_est;
                    return;
                }

                // Otherwise, choose mu base on our current strategy
                switch(cstrat) {
                case CentralityStrategy::Constant:
                    // Do a simple reduction
                    mu=sigma*mu_est;
                    break;

                case CentralityStrategy::StairStep: {

                    // Find the norm of the gradient used in the stopping
                    // criteria
                    X_Vector grad_stop(X::init(grad));
                        f_mod.grad_stop(x,grad,grad_stop);
                    Real const norm_grad=sqrt(X::innr(grad_stop,grad_stop));

                    // If we're on the first iteration, just do a simple
                    // reduction strategy
                    if(iter==1) 
                        mu=sigma*mu_est;

                    // Alternatively, if the amount of reduction in the gradient
                    // does not exceed the amount of reduction in the interior
                    // point estimate and we don't yet satisfy the gradient
                    // stopping condition, keep the interior point method at
                    // the level of the current estimate
                    else if(
                        (log10(norm_gradtyp)-log10(norm_grad)
                            < log10(mu_typ)-log10(mu_est))
                        && norm_grad >= eps_grad*norm_gradtyp
                    )
                        mu=mu_est;

                    // Otherwise, do a simple reduction
                    else
                        mu=sigma*mu_est;
                    break;

                } case CentralityStrategy::PredictorCorrector:
                    // If we're on the first iteration and never computed
                    // an objective before, do a centrality step.
                    if(f_x != f_x) 
                        mu=mu_est;

                    // Otherwise, alternate iterations between taking a
                    // centrality step and an optimality step.
                    else {
                        Real sigma0 = iter % 2 ? Real(0.) : Real(1.);
                        mu=sigma0*mu_est;
                    }
                    break;
                }
            }

           
            // Adjust the stopping conditions unless the criteria below are
            // satisfied.
            static void adjustStoppingConditions(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                Real const & mu_est=state.mu_est;
                Real const & mu_typ=state.mu_typ;
                Real const & eps_mu=state.eps_mu;
                CentralityStrategy::t const & cstrat=state.cstrat;
                Natural const & iter=state.iter;
                StoppingCondition::t & opt_stop=state.opt_stop;

                // If the estimated interior point paramter is negative, exit
                if(mu_est < Real(0.)) {
                    opt_stop=StoppingCondition::InteriorPointInstability;
                    return;
                }

                // If we're doing a predicted-corrector method, don't exit
                // on the prediction step.  It ignores the interior point
                // parameter and we really want that to be small.
                if(cstrat == CentralityStrategy::PredictorCorrector &&
                    iter % 2 == 0
                ) {
                    opt_stop=StoppingCondition::NotConverged;
                    return;
                }
                
                // Prevent convergence unless mu has been reduced to
                // eps_mu * mu_typ.
                if( opt_stop==StoppingCondition::RelativeGradientSmall &&
                    !(mu_est <= mu_typ*eps_mu) 
                ) {
                    opt_stop=StoppingCondition::NotConverged;
                    return;
                }
            }

            // Conduct a line search that preserves positivity of both the
            // primal and dual variables.  For trust-region methods, this
            // is pretty straightforward since we only need to modify our
            // step.  For a line-search algorithm, we modify the line-search
            // step length so that the farthest the line-search will look
            // is within the safe region for positivity.
            static void positivityLineSearchPrimalDual(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                Real const & gamma=state.gamma;
                AlgorithmClass::t const & algorithm_class
                    =state.algorithm_class;
                Z_Vector const & z=state.z;
                X_Vector const & x=state.x;
                Z_Vector const & h_x=state.h_x;
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                X_Vector & dx=state.dx;
                Z_Vector & dz=state.dz;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is alpha0 dx.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_(X::init(x));
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(state.alpha0,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1(X::init(x));
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1(Z::init(z));
                    h(x_tmp1,z_tmp1);

                // z_tmp1=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                Real alpha_x=Z::srch(z_tmp1,h_x);

                // Determine how far we can go in the dual variable 

                // Find the largest alpha such that
                // alpha dz + z >=0
                Real alpha_z=Z::srch(dz,z); 

                // Figure out how much to shorten the steps, if at all
                Real beta_x = alpha_x*gamma>Real(1.) ? Real(1.) : alpha_x*gamma;
                Real beta_z = alpha_z*gamma>Real(1.) ? Real(1.) : alpha_z*gamma;

                // Shorten the inequality multiplier step
                Z::scal(beta_z,dz);

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(beta_x,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= beta_x;
            }
            
            // Conduct a line search that preserves positivity of both the
            // primal and dual variables.  For trust-region methods, this
            // is pretty straightforward since we only need to modify our
            // step.  For a line-search algorithm, we modify the line-search
            // step length so that the farthest the line-search will look
            // is within the safe region for positivity.  This differs from
            // the function above since it assumes that primal and dual
            // variables are linked.
            static void positivityLineSearchPrimalDualLinked(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                Real const & gamma=state.gamma;
                Real const & mu=state.mu;
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                Z_Vector const & z=state.z;
                X_Vector const & x=state.x;
                Z_Vector const & h_x=state.h_x;
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                X_Vector & dx=state.dx;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is alpha0 dx.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_(X::init(x));
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(state.alpha0,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1(X::init(x));
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1(Z::init(z));
                h(x_tmp1,z_tmp1);

                // z_tmp2=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                Real alpha1=Z::srch(z_tmp1,h_x);

                // Determine how far we can go in the dual variable

                // z_tmp1=h'(x)dx
                h.p(x,dx_,z_tmp1);
                
                // z_tmp2 = h'(x)dx o z
                Z_Vector z_tmp2(Z::init(z));
                    Z::prod(z_tmp1,z,z_tmp2);

                // z_tmp2 = -h'(x)dx o z
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 = e
                Z::id(z_tmp1);

                // z_tmp1 = mu e
                Z::scal(mu,z_tmp1);

                // Find the largest alpha such that
                // alpha (-h'(x)dx o z) + mu e >=0
                Real alpha2=Z::srch(z_tmp2,z_tmp1);

                // Determine the farthest we can go in both variables
                Real alpha0;

                // Only the dual step is restrictive
                if( alpha1 > std::numeric_limits <Real>::max() &&
                    alpha2 <= std::numeric_limits <Real>::max() 
                )
                    alpha0 = alpha2;

                // Only the primal step is restrictive
                else if(alpha1 <= std::numeric_limits <Real>::max() &&
                    alpha2 > std::numeric_limits <Real>::max()
                )
                    alpha0 = alpha1;

                // Neither step is restrictive
                else if( alpha1 > std::numeric_limits <Real>::max() &&
                    alpha2 > std::numeric_limits <Real>::max()
                )
                    alpha0 = alpha1; 

                // Both steps are restrictive
                else
                    alpha0 = alpha1 < alpha2 ? alpha1 : alpha2;
                    
                // Next, determine if we need to back off from the
                // boundary or leave the step unchanged.
                alpha0 = alpha0*gamma > Real(1.) ?  Real(1.) : alpha0*gamma;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(alpha0,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= alpha0;

            }
            // Conduct a line search that preserves positivity of the
            // primal variable. 
            static void positivityLineSearchLogBarrier(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts 
                Real const & gamma=state.gamma;
                AlgorithmClass::t const & algorithm_class
                    =state.algorithm_class;
                X_Vector const & x=state.x;
                Z_Vector const & z=state.z;
                Z_Vector const & h_x=state.h_x;
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                X_Vector & dx=state.dx;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is alpha0 dx.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_(X::init(x));
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(state.alpha0,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1(X::init(x));
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1(Z::init(z));
                    h(x_tmp1,z_tmp1);

                // z_tmp1=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                Real alpha_x=Z::srch(z_tmp1,h_x);

                // Figure out how much to shorten the steps, if at all
                Real beta_x = alpha_x*gamma>Real(1.) ? Real(1.) : alpha_x*gamma;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(beta_x,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= beta_x;
            }
            
            // This adds the interior point through use of a state manipulator.
            template <typename ProblemClass>
            struct InteriorPointManipulator
                : public StateManipulator <ProblemClass>
            {
            private:
                // A reference to the user-defined state manipulator
                StateManipulator<ProblemClass> const & smanip;

            public:
                InteriorPointManipulator(
                    StateManipulator <ProblemClass> const & smanip_
                ) : smanip(smanip_) {}


                // Application
                void operator () (
                    typename ProblemClass::Functions::t const & fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t const & loc
                ) const {
                    // Call the user define manipulator
                    smanip(fns_,state_,loc);

                    // Dynamically cast the incoming state and fns to the
                    // to work with the interior-point spaces.  In theory,
                    // this should always work since we're doing this trickery
                    // internally.  Basically, this is required since we're
                    // inserting into either the unconstrained or equality
                    // constrained code.  Within this code, the state
                    // manipulator is hard coded to use the state for the
                    // appropriate problem even though this state is really
                    // an inequality constraint when called using the routines
                    // below.
                    typename Functions::t const & fns
                        =dynamic_cast <typename Functions::t const &> (fns_);
                    typename State::t & state 
                        =dynamic_cast <typename State::t &> (state_);

                    // Create some shorcuts
                    VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                    X_Vector const & x=state.x;
                    InteriorPointMethod::t const & ipm=state.ipm;
                    Real const & mu_est = state.mu_est;
                    Z_Vector & z=state.z;
                    Z_Vector & h_x=state.h_x;
                    Z_Vector & dz=state.dz;
                    Real & mu_typ = state.mu_typ;

                    switch(loc){
                    case OptimizationLocation::BeforeInitialFuncAndGrad:

                        // Initialize the value h(x)
                        h(x,h_x);
                
                        // Set z to be m / <h(x),e> e.  In this way,
                        // mu_est = <h(x),z> / m = 1.
                        Z::id(z);
                        Z::scal(Z::innr(z,z)/Z::innr(h_x,z),z);

                        // Estimate the interior point parameter
                        estimateInteriorPointParameter(fns,state);

                        // Set the typical value for mu
                        mu_typ=mu_est;

                        // Find an initial interior point parameter
                        findInteriorPointParameter(fns,state);

                        // In a log-barrier method, find the initial Lagrange
                        // multiplier.
                        if(ipm==InteriorPointMethod::LogBarrier)
                            findInequalityMultiplierLogBarrier(fns,state);
                        break;

                    // Adjust our step or potential step to preserve positivity
                    case OptimizationLocation::BeforeLineSearch:
                    case OptimizationLocation::BeforeActualVersusPredicted:
                        // Do the linesearch
                        switch(ipm){
                        case InteriorPointMethod::PrimalDual:
                            findInequalityMultiplierStep(fns,state);
                            positivityLineSearchPrimalDual(fns,state);
                            break;

                        case InteriorPointMethod::PrimalDualLinked:
                            positivityLineSearchPrimalDualLinked(fns,state);
                            break;
                        
                        case InteriorPointMethod::LogBarrier:
                            positivityLineSearchLogBarrier(fns,state);
                            break;
                        }
                        break;

                    // After we reject a step, make sure that we take a zero
                    // step in the inequality multiplier.  This is important
                    // in case we exit early due to small steps.
                    case OptimizationLocation::AfterRejectedTrustRegion:
                    case OptimizationLocation::AfterRejectedLineSearch:
                        Z::zero(dz);
                        break;

                    case OptimizationLocation::BeforeStep:
                        // Find the new inequality multiplier or step
                        switch(ipm){
                        case InteriorPointMethod::PrimalDual:
                            Z::axpy(Real(1.),dz,z);

                            // In theory, we start symmetric and make sure our
                            // steps are symmetric.  However, in the interest
                            // in never having a nonsymmetric dual variable,
                            // we force symmetrization here.
                            Z::symm(z);
                            break;
                        case InteriorPointMethod::PrimalDualLinked:
                            findInequalityMultiplierLinked(fns,state);
                            break;
                        case InteriorPointMethod::LogBarrier:
                            // Wait until after we update the interior
                            // point parameter to find the multiplier.
                            break;
                        }
                        break;

                    case OptimizationLocation::AfterStepBeforeGradient:
                        // Updated our cached copy of h(x)
                        h(x,h_x);

                        // Update the interior point estimate
                        estimateInteriorPointParameter(fns,state);
                        break;

                    case OptimizationLocation::AfterGradient:
                        // Update the interior point parameter
                        findInteriorPointParameter(fns,state);

                        // Update the inequality multiplier in a log-barrier
                        // method
                        if(ipm==InteriorPointMethod::LogBarrier)
                            findInequalityMultiplierLogBarrier(fns,state);
                        break; 

                    // Adjust the interior point parameter and insure that
                    // we do not converge unless the interior point parameter
                    // is small.
                    case OptimizationLocation::EndOfOptimizationIteration:
                        adjustStoppingConditions(fns,state);
                        break;

                    default:
                        break;
                    }
                }
            };

            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                StateManipulator <InequalityConstrained <Real,XX,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                StateManipulator <InequalityConstrained <Real,XX,ZZ> > const &
                    smanip,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Adds the output pieces to the state manipulator 
                DiagnosticManipulator <InequalityConstrained <Real,XX,ZZ> >
                    dmanip(smanip,msg);

                // Add the interior point pieces to the state manipulator
                InteriorPointManipulator <InequalityConstrained <Real,XX,ZZ> >
                    ipmanip(dmanip);

                // Insures that we can interact with unconstrained code
                ConversionManipulator
                    <InequalityConstrained<Real,XX,ZZ>,Unconstrained <Real,XX> >
                    cmanip(ipmanip);
                
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);
                
                // Check the inputs to the optimization
                State::check(msg,state);
                
                // Minimize the problem
                Unconstrained <Real,XX>::Algorithms
                    ::getMin_(msg,cmanip,fns,state);
            }
        };
    };
        
    // Routines that manipulate and support problems of the form
    // problem of the form
    // 
    // min_{x \in X} f(x) st g(x) = 0, h(x) >=_K 0
    //
    // where f : X -> R, g : X -> Y, and h : X -> Z
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class YY,
        template <typename> class ZZ
    > 
    struct Constrained {
        // Disallow constructors
        NO_CONSTRUCTORS(Constrained);

        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State);
            
            // Internal state of the optimization
            struct t: 
                public EqualityConstrained <Real,XX,YY>::State::t,
                public InequalityConstrained <Real,XX,ZZ>::State::t
            {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t);

                // Initialization constructors
                explicit t(
                    X_Vector const & x,Y_Vector const & y,Z_Vector const & z
                ) : Unconstrained <Real,XX>::State::t(x), 
                    EqualityConstrained <Real,XX,YY>::State::t(x,y),
                    InequalityConstrained <Real,XX,ZZ>::State::t(x,z)
                {}
            };
            
            // Check that we have a valid set of parameters.
            static void check(Messaging const & msg,t const & state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
                InequalityConstrained <Real,XX,ZZ>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {
            // Disallow constructors
            NO_CONSTRUCTORS(Restart);
       
            // Holds restart information
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <Natural> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Y_Vector> > Y_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Z_Vector> > Z_Vectors;
            
            // Checks whether we have a valid real label.
            struct is_real :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename EqualityConstrained <Real,XX,YY>::Restart
                            ::is_real()(name) ||
                        typename InequalityConstrained <Real,XX,ZZ>::Restart
                            ::is_real()(name)
                    )
                        return true;
                    else
                        return false;
                    }
            };
            
            // Checks whether we have a valid natural number label.
            struct is_nat :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename EqualityConstrained <Real,XX,YY>::Restart
                            ::is_nat()(name) ||
                        typename InequalityConstrained <Real,XX,ZZ>::Restart
                            ::is_nat()(name)
                    )
                        return true;
                    else
                        return false;
                }
            };
           
            // Checks whether we have a valid parameter label.
            struct is_param :
                public std::unary_function<std::string const &, bool>
            {
                bool operator () (std::string const & name) const {
                    if( typename EqualityConstrained <Real,XX,YY>::Restart
                            ::is_param()(name) ||
                        typename InequalityConstrained <Real,XX,ZZ>::Restart
                            ::is_param()(name)
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid variable label
            struct is_x : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( typename EqualityConstrained <Real,XX,YY>::Restart
                            ::is_x()(name) ||
                        typename InequalityConstrained <Real,XX,ZZ>::Restart
                            ::is_x()(name)
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid equality multiplier label
            struct is_y : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( typename EqualityConstrained <Real,XX,YY>::Restart
                            ::is_y()(name)
                    ) 
                        return true;
                    else
                        return false;
                }
            };
            
            // Checks whether we have a valid inequality multiplier label
            struct is_z : public std::unary_function<std::string, bool> {
                bool operator () (std::string const & name) const {
                    if( typename InequalityConstrained <Real,XX,ZZ>::Restart
                            ::is_z()(name)
                    ) 
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have valid labels
            static void checkLabels(
                Messaging const & msg,
                Reals const & reals,
                Nats const & nats,
                Params const & params,
                X_Vectors const & xs,
                Y_Vectors const & ys,
                Z_Vectors const & zs
            ) {
                Utility::checkLabels(msg,is_real(),reals.first," real name: ");
                Utility::checkLabels(msg,is_nat(),nats.first," natural name: ");
                Utility::checkLabels(
                    msg,is_param(),params.first," paramater name: ");
                Utility::checkLabels(msg,is_x(),xs.first," variable name: ");
                Utility::checkLabels(
                    msg,is_y(),ys.first," equality multiplier name: ");
                Utility::checkLabels(
                    msg,is_z(),zs.first," inequality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string const &,std::string const &,std::string>
            {
                std::string operator() (
                    std::string const & label,
                    std::string const & val
                ) const {

                    // Create a base message
                    const std::string base
                        ="During serialization, found an invalid ";

                    // Used to build the message 
                    std::stringstream ss;

                    // Check the equality parameters
                    if(typename EqualityConstrained <Real,XX,YY>
                        ::Restart::is_param()(label)
                    ) {
                        ss << typename EqualityConstrained <Real,XX,YY>
                            ::Restart::checkParamVal()(label,val);

                    // Check the inequality parameters
                    } else if (typename InequalityConstrained <Real,XX,ZZ>
                        ::Restart::is_param()(label)
                    ) {
                        ss << typename InequalityConstrained <Real,XX,ZZ>
                            ::Restart::checkParamVal()(label,val);
                    }

                    return ss.str();
                }
            };
            
            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Z_Vectors & zs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>
                    ::Restart::stateToVectors(state,xs);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::stateToVectors(state,xs,ys);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::stateToVectors(state,xs,zs);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::stateToScalars(state,reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::stateToScalars(state,reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::stateToScalars(state,reals,nats,params);
            }

            // Capture data from structures controlled by the user.  
            static void capture(
                Messaging const & msg,
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Z_Vectors & zs,
                Reals & reals,
                Nats & nats,
                Params & params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,ys,zs);

                // Check the strings used to represent parameters
                Utility::checkParams(msg,checkParamVal(),params);

                // Copy in the variables 
                Unconstrained <Real,XX>
                    ::Restart::vectorsToState(state,xs);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::vectorsToState(state,xs,ys);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::vectorsToState(state,xs,zs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>
                    ::Restart::scalarsToState(state,reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::Restart::scalarsToState(state,reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::Restart::scalarsToState(state,reals,nats,params);

                // Check that we have a valid state 
                State::check(msg,state);
            }
        };

        // All the functions required by an optimization algorithm.  Note, this
        // routine owns the memory for these operations.  
        struct Functions {
            // Disallow constructors
            NO_CONSTRUCTORS(Functions);

            // Actual storage of the functions required
            struct t: 
                public EqualityConstrained <Real,XX,YY>::Functions::t,
                public InequalityConstrained <Real,XX,ZZ>::Functions::t
            {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t);

                // Initialize all of the pointers to null
                t() : EqualityConstrained <Real,XX,YY>::Functions::t(), 
                    InequalityConstrained <Real,XX,ZZ>::Functions::t() {}
            };

            // Check that all the functions are defined
            static void check(Messaging const & msg,t const & fns) {
                EqualityConstrained <Real,XX,YY>::Functions::check(msg,fns);
                InequalityConstrained <Real,XX,ZZ>::Functions::check(msg,fns);
            }

            // Initialize any missing functions 
            static void init(
                Messaging const & msg,
                typename State::t & state,
                t& fns
            ) {
                Unconstrained <Real,XX>
                    ::Functions::init_(msg,state,fns);
                EqualityConstrained <Real,XX,YY>
                    ::Functions::init_(msg,state,fns);
                InequalityConstrained <Real,XX,ZZ>
                    ::Functions::init_(msg,state,fns);
            }
        };
        
        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
            // Disallow constructors
            NO_CONSTRUCTORS(Printer);

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { 
            }
            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>
                    ::Printer::getStateHeader_(state,out);
                EqualityConstrained <Real,XX,YY>
                    ::Printer::getStateHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>
                    ::Printer::getStateHeader_(state,out);
            }

            // Combines all of the state information
            static void getState(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                bool const & noiter,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getState_(fns,state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getState_(fns,state,blank,out);
            }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getKrylovHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylovHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylovHeader_(state,out);
            }

            // Combines all of the Krylov information
            static void getKrylov(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getKrylov_(state,blank,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylov_(state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };
        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms);

            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                StateManipulator <Constrained <Real,XX,YY,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                StateManipulator <Constrained <Real,XX,YY,ZZ> > const & smanip,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Adds the output pieces to the state manipulator 
                DiagnosticManipulator <Constrained <Real,XX,YY,ZZ> >
                    dmanip(smanip,msg);

                // Add the interior point pieces to the state manipulator
                typename InequalityConstrained <Real,XX,ZZ>::Algorithms
                    ::template InteriorPointManipulator
                    <Constrained <Real,XX,YY,ZZ> >
                    ipmanip(dmanip);

                // Add the composite step pieces to the state manipulator
                typename EqualityConstrained <Real,XX,YY>::Algorithms
                    ::template CompositeStepManipulator
                    <Constrained <Real,XX,YY,ZZ> >
                    csmanip(ipmanip,msg);

                // Insures that we can interact with unconstrained code
                ConversionManipulator
                    <Constrained<Real,XX,YY,ZZ>,Unconstrained <Real,XX> >
                    cmanip(csmanip);
                
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);
                
                // Check the inputs to the optimization
                State::check(msg,state);
                
                // Minimize the problem
                Unconstrained <Real,XX>::Algorithms
                    ::getMin_(msg,cmanip,fns,state);
            }
        };
    };
}
#endif
