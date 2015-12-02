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
#include "optizelle/linalg.h"

//---Optizelle0---
namespace Optizelle{
//---Optizelle1---

    //---ScalarValuedFunction0---
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
        virtual Real eval(Vector const & x) const = 0;

        // grad = grad f(x) 
        virtual void grad(Vector const & x,Vector & grad) const = 0;

        // H_dx = hess f(x) dx 
        virtual void hessvec(Vector const & x,Vector const & dx,Vector & H_dx)
            const = 0;

        // Allow a derived class to deallocate memory
        virtual ~ScalarValuedFunction() {}
    };
    //---ScalarValuedFunction1---
    
    template <
        typename Real,
        template <typename> class XX
    >
    struct ScalarValuedFunctionModifications {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

        // Disallow constructors
        NO_COPY_ASSIGNMENT(ScalarValuedFunctionModifications)

        // Use an empty default constructor
        ScalarValuedFunctionModifications() = default; 

        // Allow derived classes to deallocate memory
        virtual ~ScalarValuedFunctionModifications() = default; 

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


    //---VectorValuedFunction0---
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
        virtual void eval(X_Vector const & x,Y_Vector & y) const = 0;

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
    //---VectorValuedFunction1---

    //---Messaging0---
    // Defines how we output messages to the user
    struct Messaging {
        // Prints a message
        virtual void print(std::string const & msg) const;

        // Prints an error
        virtual void error(std::string const & msg) const;

        // Allow a derived class to deallocate memory
        virtual ~Messaging();
    };
    //---Messaging1---

    // A safeguard search used primarily for inequality constraints
    template <typename Real,template <typename> class XX>
    using Safeguard = std::function <
        Real(
            typename XX <Real>::Vector const & dx_base,
            typename XX <Real>::Vector const & dx_dir,
            Real const & zeta
        )>;

    // Communicates whether or not the gradient used for the step calculation
    // has been modified.  We use this for modifying the interior point
    // parameter just prior to step calculation.
    template <typename Real,template <typename> class XX>
    using GradStepModification = std::function <bool(
        typename XX <Real>::Vector const & grad_step,
        Real const & gx_reduction,
        bool const & gx_converged)>;

    // Communicates whether we need to do a second equality multiplier solve
    // for the equality constrained problem.  In an interior point method, we
    // would have modified the inequality multiplier after the globalization,
    // which means that the current multiplier is likely old and needs to be
    // recomputed.
    typedef std::function<bool()> MultiplierSolve;

    // Determines whether we're using an aboslute or relative tolerance 
    template <typename Real>
    using ToleranceSelector = std::function<Real(Real const & typ)>;

    // Which algorithm class do we use
    namespace AlgorithmClass{
        enum t : Natural{
            //---AlgorithmClass0---
            TrustRegion,            // Trust-Region algorithms
            LineSearch,             // Line-search algorithms
            UserDefined             // User provides the iterate 
            //---AlgorithmClass1---
        };

        // Converts the algorithm class to a string
        std::string to_string(t const & algorithm_class);
        
        // Converts a string to an algorithm class 
        t from_string(std::string const & algorithm_class);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // Reasons why we stop the algorithm
    namespace StoppingCondition{
        enum t : Natural{
            //---StoppingCondition0---
            NotConverged,            // Algorithm did not converge
            GradientSmall,           // Gradient was sufficiently small
            StepSmall,               // Change in the step is small
            MaxItersExceeded,        // Maximum number of iterations exceeded
            InteriorPointInstability,// Instability in the interior point method
            UserDefined              // Some user defined stopping condition 
            //---StoppingCondition1---
        };

        // Converts the stopping condition to a string 
        std::string to_string(t const & opt_stop);

        // Converts a string to a stopping condition
        t from_string(std::string const & opt_stop);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // Various operators for both Hessian approximations and preconditioners
    namespace Operators{
        enum t : Natural{
            //---Operators0---
            Identity,          // Identity approximation
            ScaledIdentity,    // Scaled identity approximation
            BFGS,              // BFGS approximation
            InvBFGS,           // Inverse BFGS approximation
            SR1,               // SR1 approximation
            InvSR1,            // Inverse SR1 approximation
            UserDefined        // User defined operator (such as the true
                               // Hessian for Newton's method)
            //---Operators1---
        };
        
        // Converts the operator type to a string 
        std::string to_string(t const & op);
        
        // Converts a string to a operator 
        t from_string(std::string const & op);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // Different kinds of search directions 
    namespace LineSearchDirection{
        enum t : Natural{
            //---LineSearchDirection0---
            SteepestDescent,          // SteepestDescent 
            FletcherReeves,           // Fletcher-Reeves CG
            PolakRibiere,             // Polak-Ribiere CG
            HestenesStiefel,          // HestenesStiefel CG
            BFGS,                     // Limited-memory BFGS 
            NewtonCG                  // Newton-CG
            //---LineSearchDirection1---
        };
        
        // Converts the line-search direction to a string 
        std::string to_string(t const & dir);
        
        // Converts a string to a line-search direction 
        t from_string(std::string const & dir);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // Different sorts of line searches
    namespace LineSearchKind{
        enum t : Natural{
            //---LineSearchKind0---
            GoldenSection,    // Golden-section search 
            BackTracking,     // BackTracking search 
            TwoPointA,        // Barzilai and Borwein's method A
            TwoPointB         // Barzilai and Borwein's method B
            //---LineSearchKind1---
        };
            
        // Converts the line-search kind to a string 
        std::string to_string(t const & kind);
        
        // Converts a string to a line-search kind 
        t from_string(std::string const & kind);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);

        // Determine whether or not the line-search checks the sufficient
        // decrease condition.
        bool is_sufficient_decrease(t const & kind);
    }
    
    // Different points in the optimization algorithm
    namespace OptimizationLocation{
        enum t : Natural{
            //---OptimizationLocation0---
            // Occurs at the start of the optimization function 
            BeginningOfOptimization,

            // Occurs before the initial function and gradient evaluation 
            BeforeInitialFuncAndGrad,

            // Occurs after the initial function and gradient evaluation 
            AfterInitialFuncAndGrad,
            
            // Occurs just before the main optimization loop 
            BeforeOptimizationLoop,
            //---OptimizationLocation1---
                       
            // Occurs at the beginning of the optimization loop
            BeginningOfOptimizationLoop,

            // Occurs just before we take the optimization step x+dx
            BeforeSaveOld,

            // Occurs just before we take the optimization step x+dx
            BeforeStep,

            // Occurs before we calculate our new step.
            BeforeGetStep,
            //---OptimizationLocation2---
            
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
            //---OptimizationLocation3---

            // Occurs after we update our quasi-Newton information. 
            AfterQuasi,
            
            // This occurs after we check our stopping condition.  This is
            // where the equality and inequality algorithms adjust the
            // stopping conditions.
            AfterCheckStop,

            // This occurs last in the optimization loop.  At this point,
            // we have already incremented our optimization iteration and
            // checked our stopping condition
            EndOfOptimizationIteration,

            // This occurs prior to the computation of the line search
            BeforeLineSearch,

            // This occurs after a rejected trust-region step
            AfterRejectedTrustRegion,
            //---OptimizationLocation4---

            // This occurs after a rejected line-search step
            AfterRejectedLineSearch,

            // This occurs prior to checking the predicted versus actual
            // reduction in a trust-region method.
            BeforeActualVersusPredicted,

            // This occurs at the end of a Krylov iteration
            EndOfKrylovIteration,

            // This occurs at the end of all optimization 
            EndOfOptimization
            //---OptimizationLocation5---
        };
            
        // Converts the optimization location to a string 
        std::string to_string(t const & loc);
        
        // Converts a string to a line-search kind 
        t from_string(std::string const & loc);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }
    
    // Different problem classes
    namespace ProblemClass{
        enum t : Natural{
            //---ProblemClass0---
            Unconstrained,         // Unconstrained optimization 
            EqualityConstrained,   // Equality constrained optimization 
            InequalityConstrained, // Inequality constrained optimization 
            Constrained            // Fully constrained optimization 
            //---ProblemClass1---
        };

        // Converts the problem class to a string
        std::string to_string(t const & problem_class);

        // Converts a string to a problem class 
        t from_string(std::string const & problem_class);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }
    
    // When and how often we compute our intrusive diagnostics 
    namespace DiagnosticScheme {
        enum t : Natural{
            //---DiagnosticScheme0---
            Never,              // Never compute our diagnostic checks 
            DiagnosticsOnly,    // No optimization.  Only diagnostics.
            EveryIteration      // Every iteration at the start of the iteration
            //---DiagnosticScheme1---
        };
        
        // Converts the diagnostic scheme to a string
        std::string to_string(t const & dscheme);
        
        // Converts a string to the diagnostic scheme 
        t from_string(std::string const & dscheme);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & dscheme);
    }
    
    // Different function diagnostics on the optimization functions 
    namespace FunctionDiagnostics {
        enum t : Natural{
            //---FunctionDiagnostics0---
            NoDiagnostics,      // No diagnostic checks
            FirstOrder,         // First-order function checks
            SecondOrder         // Second-order function checks
            //---FunctionDiagnostics1---
        };
        
        // Converts the diagnostic checks to a string
        std::string to_string(t const & diag);
        
        // Converts a string to the diagnostic checks 
        t from_string(std::string const & diag);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }
    
    // Different diagnostics on the algebra 
    namespace VectorSpaceDiagnostics {
        enum t : Natural{
            //---VectorSpaceDiagnostics0---
            NoDiagnostics,      // No diagnostic checks
            Basic,              // Test our basic vector space operations
            EuclideanJordan     // Test our Euclidean-jordan algebraic 
            //---VectorSpaceDiagnostics1---
        };
        
        // Converts the diagnostic checks to a string
        std::string to_string(t const & diag);
        
        // Converts a string to the diagnostic checks 
        t from_string(std::string const & diag);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }
    
    // Different kinds of stopping tolerances 
    namespace ToleranceKind {
        enum t : Natural{
            //---ToleranceKind0---
            Relative,           // Relative stopping tolerances
            Absolute,           // Absolute stopping tolerances 
            //---ToleranceKind1---
        };
        
        // Converts the diagnostic checks to a string
        std::string to_string(t const & eps_rel);
        
        // Converts a string to the diagnostic checks 
        t from_string(std::string const & eps_rel);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & eps_rel);
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
            Real obj_xpes=f.eval(x_op_dx);

            // f(x-eps dx)
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            Real obj_xmes=f.eval(x_op_dx);

            // f(x+2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(2.*epsilon),dx,x_op_dx);
            Real obj_xp2es=f.eval(x_op_dx);

            // f(x-2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.*epsilon),dx,x_op_dx);
            Real obj_xm2es=f.eval(x_op_dx);

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
            f.eval(x_op_dx,f_x_op_dx);
            Y::axpy(Real(8.),f_x_op_dx,dd);

            // f(x-eps dx)
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            f.eval(x_op_dx,f_x_op_dx);
            Y::axpy(Real(-8.),f_x_op_dx,dd);

            // f(x+2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(2.)*epsilon,dx,x_op_dx);
            f.eval(x_op_dx,f_x_op_dx);
            Y::axpy(Real(-1.),f_x_op_dx,dd);

            // f(x-2 eps dx)
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.)*epsilon,dx,x_op_dx);
            f.eval(x_op_dx,f_x_op_dx);
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
            typename XX <Real>::Vector const & dx,
            std::string const & name
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
            msg.print("Finite difference test on the gradient of " + name);
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
            typename XX <Real>::Vector const & dx,
            std::string const & name
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
            msg.print("Finite difference test on the Hessian of " + name);
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
            typename XX <Real>::Vector const & dxx,
            std::string const & name
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
            msg.print("Symmetry test on the Hessian of " + name);
            std::stringstream ss;
            ss<< "The absolute error between <H(x)dx,dxx> and <dx,H(x)dxx>: "
                << std::scientific << std::setprecision(16) << diff;
            msg.print(ss.str());
            
            // Return the absolute error in symmetry 
            return diff;
        }
        
        // Tests the symmetry of an operator comparing
        // <A dx,dxx> to <dx,A dxx>.
        template <
            typename Real,
            template <typename> class XX
        >
        Real operatorSymmetryCheck(
            Messaging const & msg,
            Operator <Real,XX,XX> const & A,
            typename XX <Real>::Vector const & dx,
            typename XX <Real>::Vector const & dxx,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Calculate A dx
            X_Vector A_dx(X::init(dx));
            A.eval(dx,A_dx);
            
            // Calculate A dxx
            X_Vector A_dxx(X::init(dx));
            A.eval(dxx,A_dxx);
            
            // Calculate < A dx,dxx>
            Real innr_Adx_dxx = X::innr(A_dx,dxx);
            
            // Calculate <dx,A dxx>
            Real innr_dx_Adxx = X::innr(dx,A_dxx);

            // Determine the absolute difference between the two.  This really
            // should be zero.
            Real diff=fabs(innr_Adx_dxx-innr_dx_Adxx);

            // Send a message with the result
            msg.print("Symmetry test on the operator " + name);
            std::stringstream ss;
            ss<< "The absolute error between <" << name << " dx,dxx> and <dx,"
                << name << " dxx>: " << std::scientific
                << std::setprecision(16) << diff;
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
            typename YY <Real>::Vector const & y,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            Y_Vector res(Y::init(y));

            // Calculate f'(x)dx 
            Y_Vector fp_x_dx(Y::init(y));
            f.p(x,dx,fp_x_dx);

            // Compute an ensemble of finite difference tests in a linear manner
            std::stringstream notice;
            notice << "Finite difference test on the derivative of " 
                << name;
            msg.print(notice.str());
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
            typename YY <Real>::Vector const & dy,
            std::string const & name
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
            msg.print("Adjoint test on the first derivative of " + name);
            std::stringstream ss;
            ss<<"The absolute err. between <" + name + "'(x)dx,dy> and <dx,"
                + name + "'(x)*dy>: "
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
            typename YY <Real>::Vector const & dy,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            X_Vector res(X::init(x));

            // Calculate (f''(x)dx)*dy
            X_Vector fpps_x_dx_dy(X::init(dx));
            f.pps(x,dx,dy,fpps_x_dx_dy);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the 2nd-derivative adjoint "
                "of " + name);
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
        
        // Checks the zero and innr operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real zero_innr(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create a zero vector 
            X_Vector zero(X::init(x));
            X::zero(zero);

            // Figure out it's norm
            Real norm = sqrt(X::innr(zero,zero));

            // Print out it's norm
            std::stringstream ss;
            ss << "The " << name << "::norm of zero(x) is: " << norm;
            msg.print(ss.str());

            // Return the actual norm 
            return norm;
        }
        
        // Checks the copy, axpy, and innr operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real copy_axpy_innr(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Make a copy of x 
            X_Vector xhat(X::init(x));
            X::copy(x,xhat);

            // xhat <- -0.5 x + xhat
            X::axpy(Real(-0.5),x,xhat);

            // xhat <- 0.5 x + xhat
            X::axpy(Real(0.5),x,xhat);

            // xhat <- -1.0 x + xhat
            X::axpy(Real(-1.0),x,xhat);

            // Figure out it's norm
            Real norm = sqrt(X::innr(xhat,xhat));

            // Print out it's norm
            std::stringstream ss;
            ss << "The " << name << "::norm of ((x-0.5x)+0.5x)-x is: " << norm;
            msg.print(ss.str());

            // Return the actual norm 
            return norm;
        }
        
        // Checks the copy, scal, and innr operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real copy_scal_innr(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Make a copy of x 
            X_Vector xhat(X::init(x));
            X::copy(x,xhat);

            // xhat <- 10.0 x 
            X::scal(Real(10.0),xhat);

            // Figure out || xhat || 
            Real xhat_norm = sqrt(X::innr(xhat,xhat));

            // Figure out || x ||
            Real x_norm = sqrt(X::innr(x,x));

            // Figure out 10 || x || - || xhat ||
            Real norm = Real(10.)*x_norm-xhat_norm;

            // Print out their norms
            std::stringstream ss;
            ss << "The value || 10 x || - 10 || x || in the " << name <<
                "::norm is: " << norm; 
            msg.print(ss.str());

            // Return the actual norm 
            return norm;
        }
        
        // Checks the id and prod operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real id_prod(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Find the identity element 
            X_Vector e(X::init(x));
            X::id(e);

            // xhat <- x o e
            X_Vector xhat(X::init(x));
            X::prod(x,e,xhat);

            // xhat <- x - (x o e)
            X::scal(Real(-1.),xhat);
            X::axpy(Real(1.0),x,xhat);

            // Figure out || xhat || 
            Real norm = sqrt(X::innr(xhat,xhat));

            // Print out it's norm 
            std::stringstream ss;
            ss << "The value || x - (x o e) || in the " << name <<
                "::norm is: " << norm; 
            msg.print(ss.str());

            // Return the actual norm 
            return norm;
        }
        
        // Checks the prod and linv operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real prod_linv(
            Messaging const & msg,
            typename XX <Real>::Vector const & x1,
            typename XX <Real>::Vector const & x2,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // xhat <- x1 o x2 
            X_Vector xhat(X::init(x1));
            X::prod(x1,x2,xhat);

            // xtild <- Linv(x1) xhat 
            X_Vector xtild(X::init(x1));
            X::linv(x1,xhat,xtild);

            // xtild <- -x2 + xtild
            X::axpy(Real(-1.0),x2,xtild);

            // Figure out || xtild || 
            Real norm = sqrt(X::innr(xtild,xtild));

            // Print out it's norm 
            std::stringstream ss;
            ss << "The value || x2 - linv(x1)(x1 o x2)) || in the " << name <<
                "::norm is: " << norm; 
            msg.print(ss.str());

            // Return the actual norm 
            return norm;
        }
        
        // Checks the id and srch operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real id_srch(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // e <- id 
            X_Vector e(X::init(x));
            X::id(e);

            // xhat <- -2.0 e
            X_Vector xhat(X::init(x));
            X::id(xhat);
            X::scal(Real(-2.0),xhat); 

            // alpha <- srch(-2.0 e,e) - 0.5
            Real alpha(X::srch(xhat,e));
            alpha -= Real(0.5);

            // Print out the result 
            std::stringstream ss;
            ss << "The value of " << name << "::srch(-2.0 e,e) - 0.5 is: "
                << alpha; 
            msg.print(ss.str());

            // Return the actual distance 
            return alpha;
        }

        // Define the function: f(x) = barr(x)
        template <
            typename Real,
            template <typename> class XX
        >
        struct barr : public Optizelle::ScalarValuedFunction <Real,XX> {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Evaluation of the barrier function
            Real eval(const X_Vector& x) const {
                // Create some type shortcuts
                typedef XX <Real> X;

                // Return the barrier function
                return X::barr(x);
            }

            // Gradient should be: Linv(x)e
            void grad(
                const X_Vector& x,
                X_Vector& grad
            ) const {
                // e <- id 
                X_Vector e(X::init(x));
                X::id(e);

                // grad <- linv(x)e
                X::linv(x,e,grad);
            }

            // Hessian-vector product
            void hessvec(
                const X_Vector& x,
                const X_Vector& dx,
                X_Vector& H_dx
            ) const {
                X::zero(H_dx);
            }
        };
        
        // Checks the linv, id, and barr operations
        template <
            typename Real,
            template <typename> class XX
        >
        Real linv_id_barr(
            Messaging const & msg,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & dx,
            std::string const & name
        ) {
            // Do a finite difference check on the barrier function 
            std::stringstream ss;
            ss << name << "::barr";
            return gradientCheck(msg,barr <Real,XX> (),x,dx,ss.str());
        }
        
        // Checks the innr, prod, and symm operations 
        template <
            typename Real,
            template <typename> class XX
        >
        Real innr_prod_symm(
            Messaging const & msg,
            typename XX <Real>::Vector const & dx,
            typename XX <Real>::Vector const & dxx,
            typename XX <Real>::Vector const & dxxx,
            typename XX <Real>::Vector const & dxxxx,
            std::string const & name
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // xhat <- symm(dx o dxx)
            X_Vector xhat(X::init(dx));
            X::prod(dx,dxx,xhat);
            X::symm(xhat);

            // xtild <- symm(dx o dxx) o dxxx
            X_Vector xtild(X::init(dx));
            X::prod(xhat,dxxx,xtild);

            // innr1 <- <symm(dx o dxx) o dxxx , dxxxx >
            Real innr1(X::innr(xtild,dxxxx));

            // xtild <- symm(dx o dxx) o dxxxx
            X::prod(xhat,dxxxx,xtild);

            // innr2 <- < dxxx, symm(dx o dxx) o dxxxx >
            Real innr2(X::innr(dxxx,xtild));

            // Print out the result 
            std::stringstream ss;
            ss << "The value <symm(dx o dxx) o dxxx, dxxxx> - "
            "<dxxx, symm(dx o dxx) o dxxxx> using " << name <<
                "::innr is: " << innr1-innr2; 
            msg.print(ss.str());

            // Return the actual distance 
            return innr1-innr2;
        }
    }

    //---StateManipulator0---
    // A function that has free reign to manipulate or analyze the state.
    template <typename ProblemClass>
    struct StateManipulator {
        // Disallow constructors
        NO_COPY_ASSIGNMENT(StateManipulator)

        // Give an empty default constructor
        StateManipulator() {} 

        // Application
        virtual void eval(
            typename ProblemClass::Functions::t const & fns,
            typename ProblemClass::State::t & state,
            OptimizationLocation::t const & loc
        ) const = 0; 

        // Allow the derived class to deallocate memory
        virtual ~StateManipulator() {}
    };
    //---StateManipulator1---

    // A state manipulator that does nothing
    template <typename ProblemClass>
    struct EmptyManipulator : public StateManipulator <ProblemClass> {
        // Disallow constructors
        NO_COPY_ASSIGNMENT(EmptyManipulator)

        // Give an empty default constructor
        EmptyManipulator() {}

        // Application
        void eval(
            typename ProblemClass::Functions::t const & fns,
            typename ProblemClass::State:: t& state,
            OptimizationLocation::t const & loc
        ) const {}
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
        NO_COPY_ASSIGNMENT(DiagnosticManipulator)

        // Create a reference to an existing manipulator 
        explicit DiagnosticManipulator(
            StateManipulator <ProblemClass> const & smanip_,
            Messaging const & msg_
        ) : smanip(smanip_), msg(msg_) {}

        // Application
        void eval(
            typename ProblemClass::Functions::t const & fns,
            typename ProblemClass::State:: t& state,
            OptimizationLocation::t const & loc
        ) const {

            // Call the internal manipulator 
            smanip.eval(fns,state,loc);

            // Create some shortcuts
            Natural const & msg_level=state.msg_level;
            DiagnosticScheme::t const & dscheme=state.dscheme;

            // In case we're only doing diagnostics
            if( dscheme==DiagnosticScheme::DiagnosticsOnly &&
                loc == OptimizationLocation::BeforeOptimizationLoop
            ) {
                ProblemClass::Diagnostics::checkVectorSpace(msg,fns,state);
                ProblemClass::Diagnostics::checkFunctions(msg,fns,state);
                ProblemClass::Diagnostics::checkLagrangian(msg,fns,state);
                return;
            }

            // If we're not just doing diagnostics
            switch(loc){
            case OptimizationLocation::BeforeOptimizationLoop:
                // Output the headers for the diagonstic information
                if(msg_level >= 1){
                    // Get the headers 
                    std::list <std::string> out;
                    ProblemClass::Diagnostics::getStateHeader(state,out);
                    if(false && msg_level >= 3)
                        ProblemClass::Diagnostics::getKrylovHeader(state,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                        
                    // Grab some initial diagnostic information
                    out.clear();
                    ProblemClass::Diagnostics
                        ::getState(fns,state,false,false,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                }
                break;
            
            case OptimizationLocation::BeginningOfOptimizationLoop:
                // Run our diagnostic checks
                if( dscheme==DiagnosticScheme::EveryIteration ) {
                    ProblemClass::Diagnostics::checkVectorSpace(msg,fns,state);
                    ProblemClass::Diagnostics::checkFunctions(msg,fns,state);
                    ProblemClass::Diagnostics::checkLagrangian(msg,fns,state);
                }
                break;
                
            // Output the overall state at the end of the optimization
            // iteration
            case OptimizationLocation::EndOfOptimizationIteration: 
            case OptimizationLocation::AfterRejectedTrustRegion:
            case OptimizationLocation::AfterRejectedLineSearch:
                if( msg_level >= 1){
                    // Get the diagonstic information
                    std::list <std::string> out;

                    // Don't print out the iteration information if
                    // we reject the step
                    if(    loc==OptimizationLocation
                            ::AfterRejectedTrustRegion
                        || loc==OptimizationLocation
                            ::AfterRejectedLineSearch
                    )
                        ProblemClass::Diagnostics
                            ::getState(fns,state,false,true,out);
                    else 
                        ProblemClass::Diagnostics
                            ::getState(fns,state,false,false,out);

                    // Print out blank Krylov information 
                    if(false && msg_level >=3)
                        ProblemClass::Diagnostics::getKrylov(state,true,out);

                    // Output the result
                    msg.print(std::accumulate (
                        out.begin(),out.end(),std::string()));
                }
                break;

            // Output information at the end of each Krylov iteration
            case OptimizationLocation::EndOfKrylovIteration:
                if(false && msg_level >= 3) {
                    // Get the diagonstic information
                    std::list <std::string> out;

                    // Print out blank state information
                    ProblemClass::Diagnostics::getState(
                        fns,state,true,false,out);

                    // Print out the Krylov information 
                    ProblemClass::Diagnostics::getKrylov(state,false,out);

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
        NO_COPY_ASSIGNMENT(ConversionManipulator)

        // Grab a copy of the internal manipulator
        explicit ConversionManipulator(
            const StateManipulator <Internal>& smanip_
        ) : smanip(smanip_) {}

        // Application
        void eval(
            const typename External::Functions::t& fns_,
            typename External::State::t& state_,
            OptimizationLocation::t const & loc
        ) const {
            const typename Internal::Functions::t& fns
                =dynamic_cast <typename Internal::Functions::t const &> (fns_);
            typename Internal::State::t& state 
                =dynamic_cast <typename Internal::State::t &> (state_);
            smanip.eval(fns,state,loc);
        }
    };

    // Defines the type for the restart packages
    template <typename T>
    struct RestartPackage {
        typedef std::pair <std::string,T> tuple;
        typedef std::list <tuple> t;
    };

    // A series of utiilty functions used by the routines below.
    namespace Utility {
        // Checks whether all the items are actually valids inputs.  If not, it 
        // throws an error.
        template <typename T> 
        void checkItems(
            Messaging const & msg,
            std::function <
                bool(typename RestartPackage <T>::tuple const &) > is_item,
            typename RestartPackage<T>::t const & items,
            std::string const & kind
        ) {
            // Create a base message
            const std::string base
                ="During serialization, found an invalid ";

            // Check the labels
            typename RestartPackage<T>::t::const_iterator item
                = find_if(items.begin(), items.end(),std::not1(is_item));

            if(item!=items.end()) {
                std::stringstream ss;
                ss << base << kind << item->first;
                msg.error(ss.str());
            }
        }

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
        std::string const blankSeparator = ".           ";

        // Figure out if we're at the absolute beginning of the
        // optimization.  We have to be a little saavy about this
        // since we could be on the first iteration, but in the
        // middle of a line-search or trust-region method and
        // still want to output things
        template <typename ProblemClass>
        bool is_opt_begin(
            typename ProblemClass::State::t const & state
        ) {
            // Create some shortcuts
            Natural const & iter=state.iter;
            Natural const & linesearch_iter=state.linesearch_iter;
            Natural const & rejected_trustregion=state.rejected_trustregion;
            AlgorithmClass::t const & algorithm_class=state.algorithm_class;

            return (iter==1) &&
                ((algorithm_class == AlgorithmClass::LineSearch &&
                    linesearch_iter==0) ||
                ((algorithm_class == AlgorithmClass::TrustRegion ||
                  algorithm_class == AlgorithmClass::UserDefined) &&
                    rejected_trustregion == 0));
        }
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
        NO_CONSTRUCTORS(Unconstrained)

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State)
                
            // Internal state of the optimization
            struct t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t)

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

                // Algorithm class
                AlgorithmClass::t algorithm_class;

                // Preconditioner for the Hessian
                Operators::t PH_type;

                // Hessian approximation
                Operators::t H_type;

                // Norm of a typical gradient
                Real norm_gradtyp;

                // Norm of a typical trial step
                Real norm_dxtyp;

                // Optimization variable 
                X_Vector x; 
                
                // Gradient of the objective
                X_Vector grad;
                
                // Optimization step 
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

                // Number of failed safe-guard steps before quitting the method
                Natural failed_safeguard_max;

                // Number of failed safeguard steps during the last iteration 
                Natural failed_safeguard;

                // Total number of failed safeguard steps
                Natural failed_safeguard_total;

                // Amount we truncate dx in order to maintain feasibility
                // with respect to the safeguard, which probably relates to
                // the inequailty constraint
                Real alpha_x;

                // Amount we truncate dx_n in order to maintain feasibility
                // with respect to the safeguard, which probably relates to
                // the inequailty constraint
                Real alpha_x_qn;

                // Kind of stopping tolerance
                ToleranceKind::t eps_kind;
                
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

                // Function diagnostics on f
                FunctionDiagnostics::t f_diag;

                // Function diagnostics on the Lagrangian
                FunctionDiagnostics::t L_diag;

                // Vector space diagnostics on X 
                VectorSpaceDiagnostics::t x_diag;

                // Diagnostic scheme 
                DiagnosticScheme::t dscheme;

                // Initialization constructors
                explicit t(X_Vector const & x_user) :
                    eps_grad(
                        //---eps_grad0---
                        1e-8
                        //---eps_grad1---
                    ),
                    eps_dx(
                        //---eps_dx0---
                        1e-8
                        //---eps_dx1---
                    ),
                    stored_history(
                        //---stored_history0---
                        0
                        //---stored_history1---
                    ),
                    history_reset(
                        //---history_reset0---
                        5
                        //---history_reset1---
                    ),
                    iter(
                        //---iter0---
                        1
                        //---iter1---
                    ),
                    iter_max(
                        //---iter_max0---
                        std::numeric_limits <Natural>::max() 
                        //---iter_max1---
                    ),
                    opt_stop(
                        //---opt_stop0---
                        StoppingCondition::NotConverged
                        //---opt_stop1---
                    ),
                    krylov_iter(
                        //---krylov_iter0---
                        0
                        //---krylov_iter1---
                    ),
                    krylov_iter_max(
                        //---krylov_iter_max0---
                        10
                        //---krylov_iter_max1---
                    ),
                    krylov_iter_total(
                        //---krylov_iter_total0---
                        0
                        //---krylov_iter_total1---
                    ),
                    krylov_orthog_max(
                        //---krylov_orthog_max0---
                        1
                        //---krylov_orthog_max1---
                    ),
                    krylov_stop(
                        //---krylov_stop0---
                        KrylovStop::RelativeErrorSmall
                        //---krylov_stop1---
                    ),
                    krylov_rel_err(
                        //---krylov_rel_err0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---krylov_rel_err1---
                    ),
                    eps_krylov(
                        //---eps_krylov0---
                        1e-2
                        //---eps_krylov1---
                    ),
                    algorithm_class(
                        //---algorithm_class0---
                        AlgorithmClass::TrustRegion
                        //---algorithm_class1---
                    ),
                    PH_type(
                        //---PH_type0---
                        Operators::Identity
                        //---PH_type1---
                    ),
                    H_type(
                        //---H_type0---
                        Operators::UserDefined
                        //---H_type1---
                    ),
                    norm_gradtyp(
                        //---norm_gradtyp0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---norm_gradtyp1---
                    ),
                    norm_dxtyp(
                        //---norm_dxtyp0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---norm_dxtyp1---
                    ),
                    x(X::init(x_user)),
                    grad(
                        //---grad0---
                        X::init(x_user)
                        //---grad1---
                    ),
                    dx(
                        //---dx0---
                        X::init(x_user)
                        //---dx1---
                    ),
                    x_old(
                        //---x_old0---
                        X::init(x_user)
                        //---x_old1---
                    ),
                    grad_old(
                        //---grad_old0---
                        X::init(x_user)
                        //---grad_old1---
                    ),
                    dx_old(
                        //---dx_old0---
                        X::init(x_user)
                        //---dx_old1---
                    ),
                    oldY(
                        //---oldY0---
                        // Empty
                        //---oldY1--- 
                    ),
                    oldS(
                        //---oldS0---
                        // Empty
                        //---oldS1--- 
                    ), 
                    f_x(
                        //---f_x0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---f_x1---
                    ),
                    f_xpdx(
                        //---f_xpdx0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---f_xpdx1---
                    ),
                    msg_level(
                        //---msg_level0---
                        1
                        //---msg_level1---
                    ),
                    failed_safeguard_max(
                        //---failed_safeguard_max0---
                        5 
                        //---failed_safeguard_max1---
                    ),
                    failed_safeguard(
                        //---failed_safeguard0---
                        0 
                        //---failed_safeguard1---
                    ),
                    failed_safeguard_total(
                        //---failed_safeguard_total0---
                        0 
                        //---failed_safeguard_total1---
                    ),
                    alpha_x(
                        //---alpha_x0---
                        std::numeric_limits <Real>::quiet_NaN() 
                        //---alpha_x1---
                    ),
                    alpha_x_qn(
                        //---alpha_x_qn0---
                        std::numeric_limits <Real>::quiet_NaN() 
                        //---alpha_x_qn1---
                    ),
                    eps_kind(
                        //---eps_kind0---
                        ToleranceKind::Absolute
                        //---eps_kind1---
                    ),
                    delta(
                        //---delta0---
                        1.
                        //---delta1---
                    ),
                    eta1(
                        //---eta10---
                        .1
                        //---eta11---
                    ),
                    eta2(
                        //---eta20---
                        .9
                        //---eta21---
                    ),
                    ared(
                        //---ared0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---ared1---
                    ),
                    pred(
                        //---pred0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---pred1---
                    ),
                    rejected_trustregion(
                        //---rejected_trustregion0---
                        0
                        //---rejected_trustregion1---
                    ),
                    alpha0(
                        //---alpha00---
                        1.
                        //---alpha01---
                    ),
                    alpha(
                        //---alpha0---
                        std::numeric_limits <Real>::quiet_NaN()
                        //---alpha1---
                    ),
                    c1(
                        //---c10---
                        1e-4
                        //---c11---
                    ),
                    linesearch_iter(
                        //---linesearch_iter0---
                        0
                        //---linesearch_iter1---
                    ),
                    linesearch_iter_max(
                        //---linesearch_iter_max0---
                        5
                        //---linesearch_iter_max1---
                    ),
                    linesearch_iter_total(
                        //---linesearch_iter_total0---
                        0
                        //---linesearch_iter_total1---
                    ),
                    eps_ls(
                        //---eps_ls0---
                        1e-2
                        //---eps_ls1---
                    ),
                    dir(
                        //---dir0---
                        LineSearchDirection::SteepestDescent
                        //---dir1---
                    ),
                    kind(
                        //---kind0---
                        LineSearchKind::GoldenSection
                        //---kind1---
                    ),
                    f_diag(
                        //---f_diag0---
                        FunctionDiagnostics::NoDiagnostics
                        //---f_diag1---
                    ),
                    L_diag(
                        //---L_diag0---
                        FunctionDiagnostics::NoDiagnostics
                        //---L_diag1---
                    ),
                    x_diag(
                        //---x_diag0---
                        VectorSpaceDiagnostics::NoDiagnostics
                        //---x_diag1---
                    ),
                    dscheme(
                        //---dscheme0---
                        DiagnosticScheme::Never
                        //---dscheme1---
                    )
                {
                        //---x0---
                        X::copy(x_user,x);
                        //---x1---
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
                if(!(
                    //---eps_grad_valid0---
                    state.eps_grad > Real(0.)
                    //---eps_grad_valid1---
                )) 
                    ss << "The tolerance for the gradient stopping condition "
                        "must be positive: eps_grad = " << state.eps_grad;
            
                // Check that the tolerance for the step length stopping
                // condition is positive
                else if(!(
                    //---eps_dx_valid0---
                    state.eps_dx > Real(0.)
                    //---eps_dx_valid1---
                ))
                    ss << "The tolerance for the step length stopping "
                        "condition must be positive: eps_dx = " << state.eps_dx;
                    
                    //---stored_history_valid0---
                    // Any 
                    //---stored_history_valid1---
                    
                    //---history_reset_valid0---
                    // Any 
                    //---history_reset_valid1---
        
                // Check that the current iteration is positive
                else if(!(
                    //---iter_valid0---
                    state.iter > 0
                    //---iter_valid1---
                ))
                    ss << "The current optimization iteration must be "
                        "positive: iter = " << state.iter;

                // Check that the maximum iteration is positive
                else if(!(
                    //---iter_max_valid0---
                    state.iter_max > 0
                    //---iter_max_valid1---
                ))
                    ss << "The maximum optimization iteration must be "
                        "positive: iter_max = " << state.iter_max;
                    
                    //---opt_stop_valid0---
                    // Any 
                    //---opt_stop_valid1---
                    
                    //---krylov_iter_valid0---
                    // Any 
                    //---krylov_iter_valid1---

                // Check that the maximum Krylov iteration is positive
                else if(!(
                    //---krylov_iter_max_valid0---
                    state.krylov_iter_max > 0
                    //---krylov_iter_max_valid1---
                ))
                    ss << "The maximum Krylov iteration must be "
                        "positive: krylov_iter_max = " << state.krylov_iter_max;
                    
                    //---krylov_iter_total_valid0---
                    // Any 
                    //---krylov_iter_total_valid1---

                // Check that the number of vectors we orthogonalize against
                // is at least 1.
                else if(!(
                    //---krylov_orthog_max_valid0---
                    state.krylov_orthog_max > 0
                    //---krylov_orthog_max_valid1---
                ))
                    ss << "The maximum number of vectors the Krylov method"
                    "orthogonalizes against must be positive: "
                    "krylov_orthog_max = " << state.krylov_orthog_max;
                    
                    //---krylov_stop_valid0---
                    // Any 
                    //---krylov_stop_valid1---

                    //---krylov_rel_err_valid0---
                    // Any 
                    //---krylov_rel_err_valid1---
                
                // Check that the stopping tolerance for the Krylov method is
                // positive
                else if(!(
                    //---eps_krylov_valid0---
                    state.eps_krylov > Real(0.)
                    //---eps_krylov_valid1---
                ))
                    ss << "The tolerance for the Krylov method stopping "
                        "condition must be positive: eps_krylov = "
                    << state.eps_krylov;
                    
                    //---algorithm_class_valid0---
                    // Any 
                    //---algorithm_class_valid1---
                    
                    //---PH_type_valid0---
                    // Any 
                    //---PH_type_valid1---
                    
                    //---H_type_valid0---
                    // Any 
                    //---H_type_valid1---

                // Check that the norm of a typical gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(!(
                    //---norm_gradtyp_valid0---
                    state.norm_gradtyp >= Real(0.)
                    || (state.iter==1 && state.norm_gradtyp!=state.norm_gradtyp)
                    //---norm_gradtyp_valid1---
                )) 
                    ss << "The norm of a typical gradient must be nonnegative: "
                        "norm_gradtyp = " << state.norm_gradtyp; 

                // Check that the norm of a typical trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(!(
                    //---norm_dxtyp_valid0---
                    state.norm_dxtyp >= Real(0.)
                    || (state.iter==1 && state.norm_dxtyp!=state.norm_dxtyp)
                    //---norm_dxtyp_valid1---
                )) 
                    ss << "The norm of a typical trial step must be "
                        "nonnegative: norm_dxtyp = " << state.norm_dxtyp; 
                    
                    //---x_valid0---
                    // Any 
                    //---x_valid1---
                    
                    //---grad_valid0---
                    // Any 
                    //---grad_valid1---
                    
                    //---dx_valid0---
                    // Any 
                    //---dx_valid1---
                    
                    //---x_old_valid0---
                    // Any 
                    //---x_old_valid1---
                    
                    //---grad_old_valid0---
                    // Any 
                    //---grad_old_valid1---
                    
                    //---dx_old_valid0---
                    // Any 
                    //---dx_old_valid1---
                    
                    //---oldY_valid0---
                    // Any 
                    //---oldY_valid1---
                    
                    //---oldS_valid0---
                    // Any 
                    //---oldS_valid1---

                // Check that the objective value isn't a NaN past
                // iteration 1
                else if(!(
                    //---f_x_valid0---
                    state.f_x == state.f_x || state.iter==1
                    //---f_x_valid1---
                ))
                    ss<< "The objective value must be a number: f_x = "
                        << state.f_x;

                // Check that the objective value at a trial step isn't
                // a NaN past iteration 1
                else if(!(
                    //---f_xpdx_valid0---
                    state.f_xpdx == state.f_xpdx || state.iter==1
                    //---f_xpdx_valid1---
                ))
                    ss << "The objective value at the trial step must be a "
                        "number: f_xpdx = " << state.f_xpdx;
                    
                    //---msg_level_valid0---
                    // Any 
                    //---msg_level_valid1---

                    // Any 

                // Check that the maximum number of failed safeguard points is
                // at least 1.
                else if(!(
                    //---failed_safeguard_max_valid0---
                    state.failed_safeguard_max >=1
                    //---failed_safeguard_max_valid1---
                ))
                    ss << "The maximum number of failed safeguard steps must "
                        "be positive: failed_safeguard_max = "
                        << state.failed_safeguard_max;

                    //---failed_safeguard_valid0---
                    // Any 
                    //---failed_safeguard_valid1---

                    //---failed_safeguard_total_valid0---
                    // Any 
                    //---failed_safeguard_total_valid1---

                    //---alpha_x_valid0---
                    // Any
                    //---alpha_x_valid1---
                    
                    //---alpha_x_qn_valid0---
                    // Any
                    //---alpha_x_qn_valid1---
                    
                    //---eps_kind_valid0---
                    // Any 
                    //---eps_kind_valid1---

                // Check that the trust-region radius is nonnegative 
                else if(!(
                    //---delta_valid0---
                    state.delta >= Real(0.)
                    //---delta_valid1---
                ))
                    ss<< "The trust-region radius must be nonnegative: delta = "
                        << state.delta; 

                // Check that the predicted vs. actual reduction tolerance
                // is between 0 and 1
                else if(!(
                    //---eta1_valid0---
                    state.eta1 > Real(0.) && state.eta1 < Real(1.)
                    //---eta1_valid1---
                ))
                    ss << "The tolerance for whether or not we accept a "
                        "trust-region step must be between 0 and 1: eta1 = "
                        << state.eta1;
                
                // Check that the other predicted vs. actual reduction tolerance
                // is between 0 and 1
                else if(!(
                    //---eta2_valid0---
                    state.eta2 > state.eta1 && state.eta2 < Real(1.)
                    //---eta2_valid1---
                ))
                    ss << "The tolerance for whether or not we increase the "
                        "trust-region radius must be between eta1 and 1: eta2 "
                        "= " << state.eta2;
                    
                    //---ared_valid0---
                    // Any 
                    //---ared_valid1---
                    
                    //---pred_valid0---
                    // Any 
                    //---pred_valid1---
                    
                    //---rejected_trustregion_valid0---
                    // Any 
                    //---rejected_trustregion_valid1---

                // Check that the base line-search step length is nonnegative 
                else if(!(
                    //---alpha0_valid0---
                    state.alpha0 >= Real(0.)
                    //---alpha0_valid1---
                ))
                    ss<<"The base line-search step length must be nonnegative: "
                        "alpha0 = " << state.alpha0;
                    
                    //---alpha_valid0---
                    // Any 
                    //---alpha_valid1---

                // Check that the sufficient decrease parameter lies between
                // 0 and 1.
                else if(!(
                    //---c1_valid0---
                    state.c1 > Real(0.) && state.c1 < Real(1.)
                    //---c1_valid1---
                ))
                    ss << "The sufficient decrease parameter must lie between "
                        "0 and 1: c1 = " << state.c1;
                    
                    //---linesearch_iter_valid0---
                    // Any 
                    //---linesearch_iter_valid1---

                // Check that the number of line-search iterations is positve 
                else if(!(
                    //---linesearch_iter_max_valid0---
                    state.linesearch_iter_max > 0
                    //---linesearch_iter_max_valid1---
                ))
                    ss << "The maximum number of line-search iterations must "
                        "be positive: linesearch_iter_max = "
                        << state.linesearch_iter_max;
                    
                    //---linesearch_iter_total_valid0---
                    // Any 
                    //---linesearch_iter_total_valid1---

                // Check that the stopping tolerance for the line-search
                // methods is positive
                else if(!(
                    //---eps_ls_valid0---
                    state.eps_ls > Real(0.)
                    //---eps_ls_valid1---
                )) 
                    ss << "The tolerance for the line-search stopping "
                        "condition must be positive: eps_ls = " << state.eps_ls;
                    
                    //---dir_valid0---
                    // Any 
                    //---dir_valid1---

                // If we're using a golden-section search, make sure we allow
                // at least two iterations.  In addition, if we're using the
                // Barzilai-Borwein two point Hessian approximation, make sure
                // the direction is steepest descent.
                else if(!(
                    //---kind_valid0---
                    (state.kind!=LineSearchKind::GoldenSection 
                        || state.linesearch_iter_max >= 2) &&
                    (state.kind!=LineSearchKind::TwoPointA ||
                        state.kind!=LineSearchKind::TwoPointB ||
                        state.dir==LineSearchDirection::SteepestDescent)
                    //---kind_valid1---
                )) {
                    ss << "When using a golden-section search, we require at "
                        "least 2 line-search iterations: linesearch_iter_max = "
                        << state.linesearch_iter_max << std::endl << std::endl;
                    
                    ss << "When using the Barzilai-Borwein two point Hessian "
                        "approximation line-search, the search direction must "
                        "be set to SteepestDescent: dir = "
                        << LineSearchDirection::to_string(state.dir);
                }
                    //---f_diag_valid0---
                    // Any 
                    //---f_diag_valid1---
                    
                    //---L_diag_valid0---
                    // Any 
                    //---L_diag_valid1---
                    
                    //---x_diag_valid0---
                    // Any 
                    //---x_diag_valid1---
                    
                    //---dscheme_valid0---
                    // Any 
                    //---dscheme_valid1---

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
            NO_CONSTRUCTORS(Restart)

            // Create some type shortcuts 
            typedef typename RestartPackage <Real>::t Reals;
            typedef typename RestartPackage <Natural>::t Naturals;
            typedef typename RestartPackage <std::string>::t Params;
            typedef typename RestartPackage <X_Vector>::t X_Vectors;

            // Checks whether we have a valid real 
            static bool is_real(
                typename RestartPackage <Real>::tuple const & item
            ) {
                if( item.first == "eps_grad" || 
                    item.first == "eps_dx" || 
                    item.first == "krylov_rel_err" || 
                    item.first == "eps_krylov" || 
                    item.first == "norm_gradtyp" || 
                    item.first == "norm_dxtyp" || 
                    item.first == "f_x" || 
                    item.first == "f_xpdx" ||
                    item.first == "alpha_x"  ||
                    item.first == "alpha_x_qn"||
                    item.first == "delta" || 
                    item.first == "eta1" || 
                    item.first == "eta2" || 
                    item.first == "ared" || 
                    item.first == "pred" || 
                    item.first == "alpha0" || 
                    item.first == "alpha" || 
                    item.first == "c1" || 
                    item.first == "eps_ls"
                ) 
                    return true;
                else
                    return false;
            }

            // Checks whether we have a valid natural number
            static bool is_nat(
                typename RestartPackage <Natural>::tuple const & item
            ) {
                if( item.first == "stored_history" ||
                    item.first == "history_reset" || 
                    item.first == "iter" || 
                    item.first == "iter_max" || 
                    item.first == "krylov_iter" || 
                    item.first == "krylov_iter_max" ||
                    item.first == "krylov_iter_total" || 
                    item.first == "krylov_orthog_max" ||
                    item.first == "msg_level" ||
                    item.first == "failed_safeguard_max" ||
                    item.first == "failed_safeguard" ||
                    item.first == "failed_safeguard_total" ||
                    item.first == "rejected_trustregion" || 
                    item.first == "linesearch_iter" || 
                    item.first == "linesearch_iter_max" ||
                    item.first == "linesearch_iter_total" 
                ) 
                    return true;
                else
                    return false;
            }
           
            // Checks whether we have a valid parameter 
            static bool is_param (
                typename RestartPackage <std::string>::tuple const & item
            ){
                if( (item.first=="algorithm_class" &&
                        AlgorithmClass::is_valid(item.second)) ||
                    (item.first=="opt_stop" &&
                        StoppingCondition::is_valid(item.second)) ||
                    (item.first=="krylov_stop" &&
                        KrylovStop::is_valid(item.second)) ||
                    (item.first=="H_type" &&
                        Operators::is_valid(item.second)) ||
                    (item.first=="PH_type" &&
                        Operators::is_valid(item.second)) ||
                    (item.first=="dir" &&
                        LineSearchDirection::is_valid(item.second)) ||
                    (item.first=="kind" &&
                        LineSearchKind::is_valid(item.second)) ||
                    (item.first=="f_diag" &&
                        FunctionDiagnostics::is_valid(item.second)) ||
                    (item.first=="L_diag" &&
                        FunctionDiagnostics::is_valid(item.second)) ||
                    (item.first=="x_diag" &&
                        VectorSpaceDiagnostics::is_valid(item.second)) ||
                    (item.first=="dscheme" &&
                        DiagnosticScheme::is_valid(item.second)) ||
                    (item.first=="eps_kind" &&
                        ToleranceKind::is_valid(item.second))
                )
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid variable
            static bool is_x(
                typename RestartPackage <X_Vector>::tuple const & item
            ) {
                if( item.first == "x" || 
                    item.first == "grad" || 
                    item.first == "dx" || 
                    item.first == "x_old" || 
                    item.first == "grad_old" || 
                    item.first == "dx_old" || 
                    item.first.substr(0,5)=="oldY_" || 
                    item.first.substr(0,5)=="oldS_" 
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have valid labels
            static void checkItems(
                Messaging const & msg,
                Reals const & reals,
                Naturals const & nats,
                Params const & params,
                X_Vectors const & xs
            ) {
                Utility::checkItems <Real> (
                    msg,is_real,reals,"real name: ");
                Utility::checkItems <Natural> (
                    msg,is_nat,nats,"natural name: ");
                Utility::checkItems <std::string> (
                    msg,is_param,params,"paramater: ");
                Utility::checkItems <X_Vector> (
                    msg,is_x,xs,"variable name: ");
            }

            
            // Copy out all variables.
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs
            ) {
                xs.emplace_back("x",std::move(state.x));
                xs.emplace_back("grad",std::move(state.grad));
                xs.emplace_back("dx",std::move(state.dx));
                xs.emplace_back("x_old",std::move(state.x_old));
                xs.emplace_back("grad_old",std::move(state.grad_old));
                xs.emplace_back("dx_old",std::move(state.dx_old));
                
                // Write out the quasi-Newton information with sequential names.
                // Note, we're padding the numbers with zeros.  Likely, this
                // scheme will break after 1 million vectors (6 digits).  Try
                // not to use that many.
                {Natural i=1;
                for(typename std::list<X_Vector>::iterator y=state.oldY.begin();
                    y!=state.oldY.end();
                    y++
                ){
                    std::stringstream ss;
                    ss << std::setfill('0') << std::setw(6) << i++;
                    xs.emplace_back("oldY_"+ss.str(),std::move(*y));
                }}

                // Write out the quasi-Newton information with sequential names
                {Natural i=1;
                for(typename std::list<X_Vector>::iterator s=state.oldS.begin();
                    s!=state.oldS.end();
                    s++
                ){
                    std::stringstream ss;
                    ss << std::setfill('0') << std::setw(6) << i++;
                    xs.emplace_back("oldS_"+ss.str(),std::move(*s));
                }}
            }
            
            // Copy out all non-variables.  This includes reals, naturals,
            // and parameters
            static void stateToScalars(
                typename State::t & state, 
                Reals & reals,
                Naturals & nats,
                Params & params
            ) {
                // Copy in all the real numbers
                reals.emplace_back("eps_grad",std::move(state.eps_grad));
                reals.emplace_back("eps_dx",std::move(state.eps_dx));
                reals.emplace_back("krylov_rel_err",
                    std::move(state.krylov_rel_err));
                reals.emplace_back("eps_krylov",std::move(state.eps_krylov));
                reals.emplace_back("norm_gradtyp",
                    std::move(state.norm_gradtyp));
                reals.emplace_back("norm_dxtyp",std::move(state.norm_dxtyp));
                reals.emplace_back("f_x",std::move(state.f_x));
                reals.emplace_back("f_xpdx",std::move(state.f_xpdx));
                reals.emplace_back("alpha_x",std::move(state.alpha_x));
                reals.emplace_back("alpha_x_qn",std::move(state.alpha_x_qn));
                reals.emplace_back("delta",std::move(state.delta));
                reals.emplace_back("eta1",std::move(state.eta1));
                reals.emplace_back("eta2",std::move(state.eta2));
                reals.emplace_back("ared",std::move(state.ared));
                reals.emplace_back("pred",std::move(state.pred));
                reals.emplace_back("alpha0",std::move(state.alpha0));
                reals.emplace_back("alpha",std::move(state.alpha));
                reals.emplace_back("c1",std::move(state.c1));
                reals.emplace_back("eps_ls",std::move(state.eps_ls));

                // Copy in all the natural numbers
                nats.emplace_back("stored_history",
                    std::move(state.stored_history));
                nats.emplace_back("history_reset",
                    std::move(state.history_reset));
                nats.emplace_back("iter",std::move(state.iter));
                nats.emplace_back("iter_max",std::move(state.iter_max));
                nats.emplace_back("krylov_iter",std::move(state.krylov_iter));
                nats.emplace_back("krylov_iter_max",
                    std::move(state.krylov_iter_max));
                nats.emplace_back("krylov_iter_total",
                    std::move(state.krylov_iter_total));
                nats.emplace_back("krylov_orthog_max",
                    std::move(state.krylov_orthog_max));
                nats.emplace_back("msg_level",std::move(state.msg_level));
                nats.emplace_back("failed_safeguard_max",
                    std::move(state.failed_safeguard_max));
                nats.emplace_back("failed_safeguard",
                    std::move(state.failed_safeguard));
                nats.emplace_back("failed_safeguard_total",
                    std::move(state.failed_safeguard_total));
                nats.emplace_back("rejected_trustregion",
                    std::move(state.rejected_trustregion));
                nats.emplace_back("linesearch_iter",
                    std::move(state.linesearch_iter));
                nats.emplace_back("linesearch_iter_max",
                    std::move(state.linesearch_iter_max));
                nats.emplace_back("linesearch_iter_total",
                    std::move(state.linesearch_iter_total));

                // Copy in all the parameters
                params.emplace_back("algorithm_class",
                    AlgorithmClass::to_string(state.algorithm_class));
                params.emplace_back("opt_stop",
                    StoppingCondition::to_string(state.opt_stop));
                params.emplace_back("krylov_stop",
                    KrylovStop::to_string(state.krylov_stop));
                params.emplace_back("H_type",
                    Operators::to_string(state.H_type));
                params.emplace_back("PH_type",
                    Operators::to_string(state.PH_type));
                params.emplace_back("dir",
                    LineSearchDirection::to_string(state.dir));
                params.emplace_back("kind",
                    LineSearchKind::to_string(state.kind));
                params.emplace_back("f_diag",
                    FunctionDiagnostics::to_string(state.f_diag));
                params.emplace_back("L_diag",
                    FunctionDiagnostics::to_string(state.L_diag));
                params.emplace_back("x_diag",
                    VectorSpaceDiagnostics::to_string(state.x_diag));
                params.emplace_back("dscheme",
                    DiagnosticScheme::to_string(state.dscheme));
                params.emplace_back("eps_kind",
                    ToleranceKind::to_string(
                        state.eps_kind));
            }

            // Copy in all variables.  This assumes that the quasi-Newton
            // information is being read in order.
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs
            ) {
                // Clear out oldY and oldS
                state.oldY.clear();
                state.oldS.clear();

                for(typename X_Vectors::iterator item = xs.begin();
                    item!=xs.end();
                    item++
                ){
                    if(item->first=="x") 
                        state.x = std::move(item->second);
                    else if(item->first=="grad")
                        state.grad = std::move(item->second);
                    else if(item->first=="dx")
                        state.dx = std::move(item->second);
                    else if(item->first=="x_old")
                        state.x_old = std::move(item->second);
                    else if(item->first=="grad_old")
                        state.grad_old = std::move(item->second);
                    else if(item->first=="dx_old")
                        state.dx_old = std::move(item->second);
                    else if(item->first.substr(0,5)=="oldY_")
                        state.oldY.emplace_back(std::move(item->second));
                    else if(item->first.substr(0,5)=="oldS_")
                        state.oldS.emplace_back(std::move(item->second));
                }
            }

            // Copy in all non-variables.  This includes reals, naturals,
            // and parameters
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Naturals & nats,
                Params & params
            ) {
                // Copy in any reals 
                for(typename Reals::iterator item = reals.begin();
                    item!=reals.end();
                    item++
                ){
                    if(item->first=="eps_grad")
                        state.eps_grad=std::move(item->second);
                    else if(item->first=="eps_dx")
                        state.eps_dx=std::move(item->second);
                    else if(item->first=="krylov_rel_err")
                        state.krylov_rel_err=std::move(item->second);
                    else if(item->first=="eps_krylov")
                        state.eps_krylov=std::move(item->second);
                    else if(item->first=="norm_gradtyp")
                        state.norm_gradtyp=std::move(item->second);
                    else if(item->first=="norm_dxtyp")
                        state.norm_dxtyp=std::move(item->second);
                    else if(item->first=="f_x")
                        state.f_x=std::move(item->second);
                    else if(item->first=="f_xpdx")
                        state.f_xpdx=std::move(item->second);
                    else if(item->first=="alpha_x")
                        state.alpha_x=std::move(item->second);
                    else if(item->first=="alpha_x_qn")
                        state.alpha_x_qn=std::move(item->second);
                    else if(item->first=="delta")
                        state.delta=std::move(item->second);
                    else if(item->first=="eta1")
                        state.eta1=std::move(item->second);
                    else if(item->first=="eta2")
                        state.eta2=std::move(item->second);
                    else if(item->first=="ared")
                        state.ared=std::move(item->second);
                    else if(item->first=="pred")
                        state.pred=std::move(item->second);
                    else if(item->first=="alpha0")
                        state.alpha0=std::move(item->second);
                    else if(item->first=="alpha")
                        state.alpha=std::move(item->second);
                    else if(item->first=="c1")
                        state.c1=std::move(item->second);
                    else if(item->first=="eps_ls")
                        state.eps_ls=std::move(item->second);
                }
            
                // Next, copy in any naturals
                for(typename Naturals::iterator item = nats.begin();
                    item!=nats.end();
                    item++
                ){
                    if(item->first=="stored_history")
                        state.stored_history=std::move(item->second);
                    else if(item->first=="history_reset")
                        state.history_reset=std::move(item->second);
                    else if(item->first=="iter")
                        state.iter=std::move(item->second);
                    else if(item->first=="iter_max")
                        state.iter_max=std::move(item->second);
                    else if(item->first=="krylov_iter")
                        state.krylov_iter=std::move(item->second);
                    else if(item->first=="krylov_iter_max")
                        state.krylov_iter_max=std::move(item->second);
                    else if(item->first=="krylov_iter_total")
                        state.krylov_iter_total=std::move(item->second);
                    else if(item->first=="krylov_orthog_max")
                        state.krylov_orthog_max=std::move(item->second);
                    else if(item->first=="msg_level")
                        state.msg_level=std::move(item->second);
                    else if(item->first=="failed_safeguard_max")
                        state.failed_safeguard_max=std::move(item->second);
                    else if(item->first=="failed_safeguard")
                        state.failed_safeguard=std::move(item->second);
                    else if(item->first=="failed_safeguard_total")
                        state.failed_safeguard_total=std::move(item->second);
                    else if(item->first=="rejected_trustregion")
                        state.rejected_trustregion=std::move(item->second);
                    else if(item->first=="linesearch_iter")
                        state.linesearch_iter=std::move(item->second);
                    else if(item->first=="linesearch_iter_max")
                        state.linesearch_iter_max=std::move(item->second);
                    else if(item->first=="linesearch_iter_total")
                        state.linesearch_iter_total=std::move(item->second);
                }
                    
                // Next, copy in any parameters 
                for(typename Params::iterator item = params.begin();
                    item!=params.end();
                    item++
                ){
                    if(item->first=="algorithm_class")
                        state.algorithm_class
                            = AlgorithmClass::from_string(item->second);
                    else if(item->first=="opt_stop")
                        state.opt_stop
                            = StoppingCondition::from_string(item->second);
                    else if(item->first=="krylov_stop")
                        state.krylov_stop=KrylovStop::from_string(item->second);
                    else if(item->first=="H_type")
                        state.H_type=Operators::from_string(item->second);
                    else if(item->first=="PH_type")
                        state.PH_type=Operators::from_string(item->second);
                    else if(item->first=="dir")
                        state.dir
                            = LineSearchDirection::from_string(item->second);
                    else if(item->first=="kind")
                        state.kind=LineSearchKind::from_string(item->second);
                    else if(item->first=="f_diag")
                        state.f_diag
                            = FunctionDiagnostics::from_string(item->second);
                    else if(item->first=="L_diag")
                        state.L_diag
                            = FunctionDiagnostics::from_string(item->second);
                    else if(item->first=="x_diag")
                        state.x_diag
                            = VectorSpaceDiagnostics::from_string(item->second);
                    else if(item->first=="dscheme")
                        state.dscheme
                            = DiagnosticScheme::from_string(item->second);
                    else if(item->first=="eps_kind")
                        state.eps_kind
                            = ToleranceKind::from_string(item->second);
                }
            }
            
            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Reals & reals,
                Naturals & nats,
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
                Naturals & nats,
                Params & params
            ) {

                // Check the user input 
                checkItems(msg,reals,nats,params,xs);

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
            NO_CONSTRUCTORS(Functions)

            // Actual storage of the functions required
            struct t{
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t)
                
                // Objective function
                std::unique_ptr <ScalarValuedFunction <Real,XX> > f;

                // Objective function modifications
                std::unique_ptr <ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Preconditioner for the Hessian of the objective
                std::unique_ptr <Operator <Real,XX,XX> > PH;

                // Safeguard search
                std::unique_ptr <Safeguard <Real,XX>> safeguard;

                // Gradient modification for the step calculation
                std::unique_ptr <GradStepModification <Real,XX>> gradmod;

                // Indicates whether another multiplier solve is required
                std::unique_ptr <MultiplierSolve> multsolve;

                // Determines whether we're using an absolute or relative
                // tolerance 
                std::unique_ptr <ToleranceSelector<Real>> absrel;

                // Initialize all of the pointers to null
                t() : f(nullptr), PH(nullptr) {}
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // The identity operator 
            struct Identity : public Operator <Real,XX,XX> {
                void eval(X_Vector const & dx,X_Vector & result) const{
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

                // Cached memory for the modifications on the gradient
                mutable X_Vector grad_step;

            public:
                ScaledIdentity(
                    typename Functions::t const & fns,
                    typename State::t const & state
                ) : f_mod(*(fns.f_mod)),
                    x(state.x),
                    grad(state.grad),
                    delta(state.delta),
                    grad_step(X::init(state.grad))
                {};

                void eval(X_Vector const & dx,X_Vector & result) const{
                    // Determine the norm of the gradient
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
                void eval(X_Vector const & dx, X_Vector & result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the BFGS Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work;
                    for(Natural i=0;i<oldY.size();i++)
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
                void eval(X_Vector const & dx,X_Vector & result) const {

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the SR1 Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work;
                    for(Natural i=0;i<oldY.size();i++)
                        work.emplace_back(std::move(X::init(dx)));
                    X_Vector yi_m_Bisi(X::init(dx));

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

                        // Determine yi-Bisi
                        X::copy(yi,yi_m_Bisi);
                        X::axpy(Real(-1.),Bisi,yi_m_Bisi);
                        
                        // Determine <yi-Bisi,dx>
                        Real inner_yimBisi_dx(X::innr(yi_m_Bisi,dx));

                        // Determine <yi-Bisi,si>
                        Real inner_yimBisi_si(X::innr(yi_m_Bisi,si));

                        // Determine <yi-Bisi,dx>/<y_i-Bisi,si>.
                        // Store in alpha
                        Real alpha=inner_yimBisi_dx/inner_yimBisi_si;

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

                            // Determine <yi-Bisi,sj>
                            Real inner_yimBisi_sj(X::innr(yi_m_Bisi,sj));

                            // Determine <yi-Bisj,sj> / <yi-Bisi,si>.
                            // Store in beta.
                            //
                            // CHECK THIS FORMULA
                            Real beta = inner_yimBisi_sj / inner_yimBisi_si;
                        
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
                void eval(X_Vector const & dx,X_Vector & result) const{

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
                void eval(X_Vector const & dx,X_Vector & result) const{
                    sr1.eval(dx,result);
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
                NO_DEFAULT_COPY_ASSIGNMENT(HessianAdjustedFunction)

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
                 Real eval(X_Vector const & x) const {
                    return f->eval(x);
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
                        H->eval(dx,H_dx);
                     else
                        f->hessvec(x,dx,H_dx);
                 }
            };

            // Don't do a safeguard search
            static Real noSafeguard(
                X_Vector const & dx_base,
                X_Vector const & dx_dir,
                Real const & zeta
            ) {
                return std::numeric_limits <Real>::infinity();
            }

            // Don't modify the gradient used for the step
            static bool noGradStepModification(
                X_Vector const & grad_step,
                Real const & gx_reduction,
                bool const & gx_converged
            ) {
                return false;
            }

            // Don't do another multiplier solve 
            static bool noMultiplierSolve() {
                return false;
            }

            // Flips between absolute and relative tolerances 
            static Real absrelSwitch(
                ToleranceKind::t const & eps_kind,
                Real const & typ
            ) {
                return eps_kind==ToleranceKind::Relative ? typ : Real(1.);
            }

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

                // Check that all functions are defined (namely, the 
                // objective).
                check(msg,fns);

                // Modify the objective function if necessary
                fns.f.reset(new HessianAdjustedFunction(msg,state,fns));

                // We don't need to safeguard our steps nor modify our gradient
                // for the step modification
                fns.safeguard = std::make_unique <Safeguard <Real,XX>>(
                    noSafeguard);
                fns.gradmod= std::make_unique <GradStepModification <Real,XX>>(
                    noGradStepModification);
                fns.absrel = std::make_unique <ToleranceSelector<Real>>(
                    std::bind(
                        absrelSwitch,
                        std::cref(state.eps_kind),
                        std::placeholders::_1));

                // No additional multiplier solves required
                fns.multsolve = std::make_unique <MultiplierSolve> (
                    noMultiplierSolve);
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
        struct Diagnostics {
            // Disallow constructors
            NO_CONSTRUCTORS(Diagnostics)

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                LineSearchDirection::t const & dir=state.dir;
                Natural const & msg_level = state.msg_level;

                // Basic information
                out.emplace_back(Utility::atos("iter"));
                out.emplace_back(Utility::atos("f(x)"));
                out.emplace_back(Utility::atos("||grad||"));
                out.emplace_back(Utility::atos("||dx||"));
                
                // More detailed information
                if(msg_level >= 2) {
                    out.emplace_back(Utility::atos("merit(x)"));

                    // In case we're using a Krylov method
                    if(    algorithm_class==AlgorithmClass::TrustRegion
                        || dir==LineSearchDirection::NewtonCG
                    ){
                        out.emplace_back(Utility::atos("kry_iter"));
                        out.emplace_back(Utility::atos("kry_err"));
                        out.emplace_back(Utility::atos("kry_stop"));
                    }

                    // In case we're using a line-search method
                    if(algorithm_class==AlgorithmClass::LineSearch) {
                        out.emplace_back(Utility::atos("ls_iter"));
                        out.emplace_back(Utility::atos("alpha0"));
                        out.emplace_back(Utility::atos("alpha"));
                    }

                    // In case we're using a trust-region method 
                    if(algorithm_class==AlgorithmClass::TrustRegion) {
                        out.emplace_back(Utility::atos("ared"));
                        out.emplace_back(Utility::atos("pred"));
                        out.emplace_back(Utility::atos("ared/pred"));
                        out.emplace_back(Utility::atos("delta"));
                    }
                }

                // Even more detailed information
                if(msg_level >= 3) {
                    // In case we're using a Krylov method
                    if(    algorithm_class==AlgorithmClass::TrustRegion
                        || dir==LineSearchDirection::NewtonCG
                    ){
                        out.emplace_back(Utility::atos("kry_itr_tot"));
                    }
                }
            }

            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getStateHeader_(
                    state,out);
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
                Natural const & krylov_iter_total=state.krylov_iter_total;
                Real const & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Natural const & linesearch_iter=state.linesearch_iter;
                Real const & alpha0=state.alpha0;
                Real const & alpha=state.alpha;
                Real const & ared=state.ared;
                Real const & pred=state.pred;
                Real const & delta=state.delta;
                AlgorithmClass::t const & algorithm_class=state.algorithm_class;
                LineSearchDirection::t const & dir=state.dir;
                Natural const & rejected_trustregion=state.rejected_trustregion;
                Natural const & msg_level=state.msg_level;

                // Figure out if we're at the absolute beginning of the
                // optimization.
                bool opt_begin = Utility::is_opt_begin <Unconstrained> (state);

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
                out.emplace_back(Utility::atos(norm_grad));
                if(!opt_begin) 
                    out.emplace_back(Utility::atos(norm_dx));
                else
                    out.emplace_back(Utility::blankSeparator);
                
                // More detailed information 
                if(msg_level >=2) {
                    out.emplace_back(Utility::atos(merit_x));

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
                        if(!opt_begin) {
                            out.emplace_back(Utility::atos(linesearch_iter));
                            out.emplace_back(Utility::atos(alpha0));
                            out.emplace_back(Utility::atos(alpha));
                        } else 
                            for(Natural i=0;i<3;i++)
                                out.emplace_back(Utility::blankSeparator);
                    }
                    
                    // In case we're using a trust-region method
                    if(algorithm_class==AlgorithmClass::TrustRegion) {
                        if(!opt_begin) {
                            out.emplace_back(Utility::atos(ared));
                            out.emplace_back(Utility::atos(pred));
                            out.emplace_back(Utility::atos(ared/pred));
                            out.emplace_back(Utility::atos(delta));
                        } else  
                            for(Natural i=0;i<4;i++)
                                out.emplace_back(Utility::blankSeparator);
                    }
                }
                
                // Even more detail 
                if(msg_level >=3) {
                    // In case we're using a Krylov method
                    if(    algorithm_class==AlgorithmClass::TrustRegion
                        || dir==LineSearchDirection::NewtonCG
                    ){
                        if(!opt_begin) {
                            out.emplace_back(Utility::atos(krylov_iter_total));
                        } else 
                            for(Natural i=0;i<1;i++)
                                out.emplace_back(Utility::blankSeparator);
                    }
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
                Unconstrained <Real,XX>::Diagnostics
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
                Unconstrained <Real,XX>::Diagnostics::getKrylovHeader_(
                    state,out);
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
                Unconstrained <Real,XX>::Diagnostics::getKrylov_(
                    state,blank,out);
            }

            // Runs the specified function diagnostics 
            static void checkFunctions_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                X_Vector const & x=state.x;
                FunctionDiagnostics::t const & f_diag=state.f_diag;
               
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 
                X_Vector dxx(X::init(x));
                    X::rand(dxx);

                // Run the diagnostics
                switch(f_diag) {
                    case FunctionDiagnostics::FirstOrder:
                        msg.print("Diagnostics on the function f");
                        Optizelle::Diagnostics::gradientCheck(msg,f,x,dx,"f");
                        msg.print("");
                        break;
                    case FunctionDiagnostics::SecondOrder:
                        msg.print("Diagnostics on the function f");
                        Optizelle::Diagnostics::gradientCheck(msg,f,x,dx,"f");
                        Optizelle::Diagnostics::hessianCheck(msg,f,x,dx,"f");
                        Optizelle::Diagnostics::hessianSymmetryCheck(
                            msg,f,x,dx,dxx,"f");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified function diagnostics 
            static void checkFunctions(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkFunctions_(
                    msg,fns,state);
            }
            
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) { }
            
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkLagrangian_(
                    msg,fns,state);
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                VectorSpaceDiagnostics::t const & x_diag=state.x_diag;
                X_Vector const & x=state.x;
               
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 

                // Run the diagnostics
                switch(x_diag) {
                    case VectorSpaceDiagnostics::Basic:
                        msg.print("Diagnostics on the vector-space X");
                        Optizelle::Diagnostics::zero_innr <Real,XX> (msg,x,"X");
                        Optizelle::Diagnostics::copy_axpy_innr <Real,XX> (
                            msg,dx,"X");
                        Optizelle::Diagnostics::copy_scal_innr <Real,XX> (
                            msg,dx,"X");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms)

            // Checks a set of stopping conditions
            static StoppingCondition::t checkStop(
                typename Functions::t const & fns, 
                typename State::t const & state
            ){
                // Create some shortcuts
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                auto const & absrel = *(fns.absrel);
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
                if(norm_dx< eps_dx * absrel(norm_dxtyp))
                    return StoppingCondition::StepSmall;
                
                // Check whether the norm is small relative to some typical
                // gradient
                if(norm_grad < eps_grad * absrel(norm_gradtyp))
                    return StoppingCondition::GradientSmall;

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

                // Allocate memory for temporaries
                mutable X_Vector H_dx;

            public:
                // Take in the objective and the base point during construction 
                HessianOperator(
                    typename Functions::t const & fns,
                    X_Vector const & x_)
                : f(*(fns.f)), f_mod(*(fns.f_mod)), x(x_), H_dx(X::init(x_))
                {}

                // Basic application
                void eval(X_Vector const & dx,X_Vector & result)
                    const
                {
                    // H_dx <- H dx
                    f.hessvec(x,dx,H_dx);

                    // result <-  (H + f_mod) dx
                    f_mod.hessvec_step(x,dx,H_dx,result);
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
                // m(dx) = f(x) + < grad,dx > + 0.5 < H(x)dx,dx >
                Real model_dx= merit_x + X::innr(grad_step,dx)
                    + Real(.5)*X::innr(Hdx_step,dx);

                // Determine x+dx
                X_Vector x_p_dx(X::init(x));
                X::copy(dx,x_p_dx);
                X::axpy(Real(1.),x,x_p_dx);

                // Determine the merit function evaluated at x+dx
                f_xpdx=f.eval(x_p_dx);
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
                auto const & gradmod = *(fns.gradmod);
                auto const & absrel = *(fns.absrel);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                auto const & safeguard = *(fns.safeguard); 
                Real const & eps_dx=state.eps_dx;
                Real const & eps_krylov=state.eps_krylov;
                Natural const & krylov_iter_max=state.krylov_iter_max;
                Natural const & krylov_orthog_max=state.krylov_orthog_max;
                Real const & delta=state.delta;
                X_Vector const & x=state.x;
                X_Vector const & grad=state.grad;
                Real const & norm_dxtyp=state.norm_dxtyp;
                auto const & failed_safeguard_max = state.failed_safeguard_max;
                Natural & rejected_trustregion=state.rejected_trustregion;
                X_Vector & dx=state.dx;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                Real & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                Natural & history_reset=state.history_reset;
                Real & alpha = state.alpha;
                Real & alpha0 = state.alpha0;
                auto & failed_safeguard = state.failed_safeguard;
                auto & failed_safeguard_total = state.failed_safeguard_total;
                auto & alpha_x = state.alpha_x;

                // Allocate a little bit of work space
                X_Vector x_tmp1(X::init(x));

                // Allocate memory for the Cauchy point
                X_Vector dx_cp(X::init(x));

                // Find -grad f(x)
                X_Vector grad_step(X::init(x));
                    f_mod.grad_step(x,grad,grad_step);
                if(gradmod(grad_step,Real(0.),true))
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector minus_grad(X::init(x));
                    X::copy(grad_step,minus_grad);
                    X::scal(Real(-1.),minus_grad);

                // Create the Hessian operator
                HessianOperator H(fns,x);

                // Continue to look for a step until one comes back as valid
                for(rejected_trustregion=0;
                    true; 
                ) {
                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::BeforeGetStep);

                    // Set the trust-region center
                    X::zero(x_tmp1);

                    // Keep track of the residual errors
                    Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                    Real residual_err(std::numeric_limits <Real>::quiet_NaN());

                    // Create the simplified safeguard function
                    auto simplified_safeguard =
                        SafeguardSimplified <Real,XX>(std::bind(
                            safeguard,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            Real(1.)));

                    // Find the trial step 
                    truncated_cg(
                        H,
                        minus_grad,
                        PH,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        delta,
                        x_tmp1,
                        false,
                        failed_safeguard_max,
                        simplified_safeguard,
                        dx,
                        dx_cp,
                        residual_err0,
                        residual_err,
                        krylov_iter,
                        krylov_stop,
                        failed_safeguard,
                        alpha_x);

                    // Calculate the Krylov error
                    krylov_rel_err = residual_err / residual_err0;
                    krylov_iter_total += krylov_iter;

                    // Keep track of the number of failed safeguard steps
                    failed_safeguard_total+=failed_safeguard;

                    // Manipulate the state if required
                    smanip.eval(fns,state,
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
                    smanip.eval(fns,state,
                        OptimizationLocation::AfterRejectedTrustRegion);

                    // Alternatively, check if the step becomes so small
                    // that we're not making progress or if we have a step
                    // with Nans in it.  In this case, break and allow the
                    // stopping conditions to terminate optimization.  We use a
                    // zero length step so that we do not modify the current
                    // iterate.
                    Real norm_dx = sqrt(X::innr(dx,dx));
                    if(norm_dx<eps_dx*absrel(norm_dxtyp) || norm_dx!=norm_dx) {
                        X::zero(dx);
                        break;
                    }
                } 
                
                // Set line-search parameters in such a way that they are
                // consistent to what just happened in the trust-region
                // method.  This helps keep this consistent if we ever switch
                // to a line-search method.
                alpha = Real(1.);
                alpha0 = delta/Real(2.);
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
                PH.eval(grad_step,dx);
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
                    PH.eval(grad_step,dx);
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
                    PH.eval(grad_step,PH_grad_step);
                X_Vector PH_grad_old_step(X::init(grad_old_step));
                    PH.eval(grad_old_step,PH_grad_old_step);

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
                    PH.eval(grad_step,PH_grad_step);
                X_Vector PH_grad_old_step(X::init(grad_old_step));
                    PH.eval(grad_old_step,PH_grad_old_step);
                    
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
                    PH.eval(grad_step,PH_grad_step);
                    
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
                Hinv.eval(grad_step,dx);

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
                Real f_mu=f.eval(x_p_dx);
                Real merit_mu=f_mod.merit(x_p_dx,f_mu);

                // lambda
                X::copy(x,x_p_dx);
                X::axpy(lambda,dx,x_p_dx);
                Real f_lambda=f.eval(x_p_dx);
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
                        f_lambda=f_mu;
                        mu=a+beta*(b-a);

                        X::copy(x,x_p_dx);
                        X::axpy(mu,dx,x_p_dx);
                        f_mu=f.eval(x_p_dx);
                        merit_mu=f_mod.merit(x_p_dx,f_mu);

                    // Otherwise, the objective is greater on the right, so
                    // bracket on the left
                    } else {
                        b=mu;
                        mu=lambda;
                        merit_mu=merit_lambda;
                        f_mu=f_lambda;
                        mu=a+beta*(b-a);
                        lambda=a+(1-beta)*(b-a);
                
                        X::copy(x,x_p_dx);
                        X::axpy(lambda,dx,x_p_dx);
                        f_lambda=f.eval(x_p_dx);
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
                f_xpdx=f.eval(x_p_adx);

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
                f_xpdx=f.eval(x_p_adx);

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
                auto const & gradmod = *(fns.gradmod);
                auto const & absrel = *(fns.absrel);
                Operator <Real,XX,XX> const & PH=*(fns.PH);
                auto const & safeguard = *(fns.safeguard); 
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
                Real const & c1=state.c1;
                auto const & failed_safeguard_max = state.failed_safeguard_max;
                X_Vector & dx=state.dx;
                Real & f_xpdx=state.f_xpdx;
                Real & alpha0=state.alpha0;
                Real & alpha=state.alpha;
                Real & krylov_rel_err=state.krylov_rel_err;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                Real & delta=state.delta;
                auto & failed_safeguard = state.failed_safeguard;
                auto & failed_safeguard_total = state.failed_safeguard_total;
                auto & alpha_x = state.alpha_x;
                
                // Manipulate the state if required
                smanip.eval(fns,state,OptimizationLocation::BeforeGetStep);
                
                // Modify the gradient if need be 
                auto grad_step = X::init(x);
                f_mod.grad_step(x,grad,grad_step);
                if(gradmod(grad_step,Real(0.),true))
                    f_mod.grad_step(x,grad,grad_step);

                // Create the trust-region offset 
                auto x_offset = X::init(x);
                X::zero(x_offset);

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
                    // Allocate memory for the Cauchy point
                    X_Vector dx_cp(X::init(x));

                    // Create the Hessian operator
                    HessianOperator H(fns,x);

                    // Find -grad f(x)
                    X_Vector grad_step(X::init(grad));
                        f_mod.grad_step(x,grad,grad_step);
                    X_Vector minus_grad(X::init(x));
                        X::copy(grad_step,minus_grad);
                        X::scal(Real(-1.),minus_grad);

                    // Keep track of the residual errors
                    Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                    Real residual_err(std::numeric_limits <Real>::quiet_NaN());

                    // Turn off the safeguards for the line-search Newton-CG 
                    auto simplified_safeguard =
                        SafeguardSimplified <Real,XX>(std::bind(
                            Safeguard <Real,XX> (Functions::noSafeguard),
                            std::placeholders::_1,
                            std::placeholders::_2,
                            Real(1.)));

                    // Find the trial step 
                    truncated_cg(
                        H,
                        minus_grad,
                        PH,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        std::numeric_limits <Real>::infinity(),
                        x_offset,
                        false,
                        failed_safeguard_max,
                        simplified_safeguard,
                        dx,
                        dx_cp,
                        residual_err0,
                        residual_err,
                        krylov_iter,
                        krylov_stop,
                        failed_safeguard,
                        alpha_x);

                    // Calculate the Krylov error
                    krylov_rel_err = residual_err / residual_err0;
                    krylov_iter_total += krylov_iter;

                    // Keep track of the number of failed safeguard steps
                    failed_safeguard_total+=failed_safeguard;
                    break;
                }}

                // Manipulate the state if required
                smanip.eval(fns,state,OptimizationLocation::BeforeLineSearch);

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

                    // Do a safeguard to insure that we don't violate some kind
                    // of constraint if one exists.  Basically, we know that we
                    // can search up to alpha0 out in front of the direction,
                    // so we make sure that x + alpha0 dx is safe.  If not, we
                    // move back alpha.
                    auto zero = X::init(x);
                    X::zero(zero);
                    auto alpha_dx = X::init(x);
                    X::copy(dx,alpha_dx);
                    X::scal(alpha0,alpha_dx);
                    alpha_x = std::min(
                        Real(1.),safeguard(zero,alpha_dx,Real(1.)));
                    if(alpha_x < Real(1.))
                        alpha0 *= alpha_x;

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
                        auto merit_pred =merit_x+c1*alpha*X::innr(grad_step,dx);
                        sufficient_decrease = 
                            (merit_xpdx==merit_xpdx) &&
                            (merit_xpdx <= merit_pred);

                        // If we've not satisfied the sufficient decrease
                        // condition, cut the step
                        if( !sufficient_decrease ) {

                            // Decrease the size of the base line-search
                            // parameter 
                            alpha0/=Real(2.);

                            // Scale the search direction for the the state
                            // manipulator, which produces output
                            X::scal(alpha,dx);

                            // Check if the step becomes so small that we're not
                            // making progress or if we have a step with NaNs in
                            // it.  In this case, take a zero step and allow
                            // the stopping conditions to exit
                            Real norm_dx = sqrt(X::innr(dx,dx));
                            if( norm_dx < eps_dx * absrel(norm_dxtyp) ||
                                norm_dx != norm_dx
                            ) {
                                X::zero(dx);
                                break;
                            }

                            // Manipulate the state if required
                            smanip.eval(fns,state,
                                OptimizationLocation::AfterRejectedLineSearch);

                            // Rescale the search direction in case we're not
                            // quite done searching yet.
                            X::scal(Real(1.0)/alpha,dx);
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
                
                // Set trust-region parameters in such a way that they are
                // consistent to what just happened in the line-search method 
                // method.  This helps keep this consistent if we ever switch
                // to a trust-region method.
                delta = Real(2.)*alpha0;
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
                    smanip.eval(fns,state,OptimizationLocation::GetStep);
                    break;
                }
            }

            // Updates the quasi-Newton information
            static void updateQuasi(
                Messaging const & msg,
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

                // If we're using SR1, check that
                // 
                // | s'(y-Bs) | > epsilson max_i | si'(yi-Bsi) |
                //
                // and
                //
                // | s'(y-Bs) | > epsilon ||s|| ||y - Bs||
                //
                // where we choose epsilon to be the square root of machine
                // precision.
                if( PH_type==Operators::InvSR1 ||
                    H_type==Operators::SR1
                ) {
                    // Bs <- B s
                    X_Vector Bs(X::init(x));
                        typename Functions::SR1(msg,state).eval(s,Bs);

                    // y_m_Bs <- y-Bs
                    X_Vector y_m_Bs(X::init(x));
                        X::copy(y,y_m_Bs);
                        X::axpy(Real(-1.),Bs,y_m_Bs);

                    // norm_s_2 = || s ||^2
                    Real norm_s_2(X::innr(s,s));

                    // norm_ymBs_2 = || y - Bs ||^2
                    Real norm_ymBs_2(X::innr(y_m_Bs,y_m_Bs));

                    // Compute a measure of how much interesting new information
                    // we'll add to the SR1 operator
                    Real innr_s_ymBs(fabs(X::innr(s,y_m_Bs)));

                    // Repeat the above step for the existing vectors
                    std::list <X_Vector> oldSS;
                    std::list <X_Vector> oldYY;
                    Real innr_si_ymBsi(0.);
                    Natural m=oldS.size();
                    for(int i=0;i<m;i++) {
                        // Remove the last vector in quasi-Newton information
                        oldSS.splice(oldSS.end(),oldS,oldS.begin());
                        oldYY.splice(oldYY.end(),oldY,oldY.begin());

                        // Bs <- B si
                        typename Functions::SR1(msg,state).eval(
                            oldSS.back(),Bs);

                        // y_m_Bs
                        X::copy(oldYY.back(),y_m_Bs);
                        X::axpy(Real(-1.),Bs,y_m_Bs);
                    
                        // Compute a measure of how much interesting new
                        // information we've already added
                        Real tmp(fabs(X::innr(oldSS.back(),y_m_Bs)));
                        innr_si_ymBsi =
                            tmp > innr_si_ymBsi ? tmp : innr_si_ymBsi;
                    }

                    // Put all the vectors back.  I'm sure there's a better
                    // way to cache this information.
                    oldS.splice(oldS.begin(),oldSS);
                    oldY.splice(oldY.begin(),oldYY);

                    // If the new vector doesn't add much, ignore it
                    if( innr_s_ymBs <=
                            sqrt(std::numeric_limits <Real>::epsilon())
                                *innr_si_ymBsi || 
                        innr_s_ymBs <= sqrt( 
                            std::numeric_limits <Real>::epsilon() *
                            norm_s_2 *
                            norm_ymBs_2)
                    )
                        return;
                }

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
                DiagnosticScheme::t const & dscheme=state.dscheme;
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
                smanip.eval(fns,state,
                    OptimizationLocation::BeginningOfOptimization);

                // Evaluate the objective function and gradient if we've not
                // done so already.  Some diagnostics may require this
                // information.
                if(f_x != f_x) {

                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation
                        ::BeforeInitialFuncAndGrad);

                    // Sometimes, we can calculate the gradient and objective
                    // simultaneously.  Hence, it's best to calculate the
                    // gradient first and then possibly cache the objective
                    f.grad(x,grad);
                    f_x=f.eval(x);
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
                    smanip.eval(fns,state,
                        OptimizationLocation::AfterInitialFuncAndGrad);
                }

                // Manipulate the state if required
                smanip.eval(fns,state,
                    OptimizationLocation::BeforeOptimizationLoop);

                // Primary optimization loop
                while(opt_stop==StoppingCondition::NotConverged &&
                      dscheme!=DiagnosticScheme::DiagnosticsOnly
                ) {

                    // Manipulate the state if required
                    smanip.eval(fns,state,
                        OptimizationLocation::BeginningOfOptimizationLoop);

                    // Get a new optimization iterate.  
                    getStep(msg,smanip,fns,state);

                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::BeforeSaveOld);

                    // Save the old variable, gradient, and trial step.  This
                    // is useful for both CG and quasi-Newton methods.
                    X::copy(x,x_old);
                    X::copy(grad,grad_old);
                    X::copy(dx,dx_old);

                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::BeforeStep);

                    // Move to the new iterate
                    X::axpy(Real(1.),dx,x);

                    // Save the size of the first step
                    if(iter==1)
                        norm_dxtyp=std::sqrt(X::innr(dx,dx));

                    // Manipulate the state if required
                    smanip.eval(fns,state,
                        OptimizationLocation::AfterStepBeforeGradient);

                    // Find the new objective value and gradient
                    f_x=f_xpdx;
                    f.grad(x,grad);
                    
                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::AfterGradient);
                    
                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::BeforeQuasi);

                    // Update the quasi-Newton information
                    updateQuasi(msg,fns,state);
                    
                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::AfterQuasi);

                    // Increase the iteration
                    iter++;
                    
                    // Check the stopping condition
                    opt_stop=checkStop(fns,state);

                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::AfterCheckStop);

                    // Manipulate the state if required
                    smanip.eval(fns,state,
                        OptimizationLocation::EndOfOptimizationIteration);
                } 
                        
                // Manipulate the state one final time if required
                smanip.eval(fns,state,OptimizationLocation::EndOfOptimization);
            }
            
            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                EmptyManipulator <Unconstrained <Real,XX> > smanip;

                // Minimize the problem
                getMin(msg,fns,state,smanip);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state,
                StateManipulator <Unconstrained <Real,XX> > const & smanip
            ){
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);

                // Check the inputs to the optimization
                State::check(msg,state);

                // Add the output to the state manipulator
                DiagnosticManipulator <Unconstrained<Real,XX> >
                    dmanip(smanip,msg);

                // Minimize the problem
                getMin_(msg,dmanip,fns,state);
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
        NO_CONSTRUCTORS(EqualityConstrained)

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
            NO_CONSTRUCTORS(State)

            // Internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t)

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

                // Number of iterations taken by the augmented system solve
                Natural augsys_qn_iter;
                Natural augsys_pg_iter;
                Natural augsys_proj_iter;
                Natural augsys_tang_iter;
                Natural augsys_lmh_iter;

                // Total number of iterations taken by the augmented system
                // solve
                Natural augsys_qn_iter_total;
                Natural augsys_pg_iter_total;
                Natural augsys_proj_iter_total;
                Natural augsys_tang_iter_total;
                Natural augsys_lmh_iter_total;
                Natural augsys_iter_total;

                // Error in the augmented system solve 
                Real augsys_qn_err;
                Real augsys_pg_err;
                Real augsys_proj_err;
                Real augsys_tang_err;
                Real augsys_lmh_err;

                // Target error in the augmented system solve 
                Real augsys_qn_err_target;
                Real augsys_pg_err_target;
                Real augsys_proj_err_target;
                Real augsys_tang_err_target;
                Real augsys_lmh_err_target;
                
                // Equality constraint evaluated at x.  We use this in the
                // quasinormal step as well as in the computation of the
                // linear Taylor series at x in the direciton dx_n.
                Y_Vector g_x;

                // A typical norm for norm_gx.  Generally, we just take
                // the value at the first iteration.
                Real norm_gxtyp;
                
                // Linear Taylor series at x in the direction dx_n.  We use  
                // this both in the predicted reduction as well as the
                // residual predicted reduction. 
                Y_Vector gpxdxn_p_gx;

                // Derivative of the constraint applied to the tangential step
                // this is used in the residual predicted reduction.
                Y_Vector gpxdxt;

                // Norm of gpxdxn_p_gx.  We use this in the penalty parameter
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
                
                // Hessian applied to the normal step.  We require this in 
                // W_gradpHdxn as well as the predicted reduction.
                X_Vector H_dxn;

                // Quantity grad f(x) + g'(x)*y + H dx_n projected into the
                // null-space of the constraints.  We require this in the
                // tangential subproblem and the predicted reduction.
                X_Vector W_gradpHdxn;
                
                // Hessian applied to the uncorrected tangential step.  We 
                // require this in the predicted reduction.
                X_Vector H_dxtuncorrected;

                // Function diagnostics on g
                FunctionDiagnostics::t g_diag;

                // Vector space diagnostics on Y 
                VectorSpaceDiagnostics::t y_diag;

                // Initialization constructors
                explicit t(X_Vector const & x_user,Y_Vector const & y_user) : 
                    Unconstrained <Real,XX>::State::t(x_user),
                    y(Y::init(y_user)),
                        //---y0---
                        // Equality constrained
                        // argmin_y || grad f(x) + g'(x)*y ||
                        //
                        // Constrained
                        // argmin_y || grad f(x) + g'(x)*y - h'(x)*z ||
                        //---y1---
                    dy(
                        //---dy0---
                        Y::init(y_user)
                        //---dy1---
                    ),
                    zeta(
                        //---zeta0---
                        0.8
                        //---zeta1---
                    ),
                    eta0(
                        //---eta00---
                        0.5
                        //---eta01---
                    ),
                    rho(
                        //---rho0---
                        1.0
                        //---rho1---
                    ),
                    rho_old(
                        //---rho_old0---
                        rho
                        //---rho_old1---
                    ),
                    rho_bar(
                        //---rho_bar0---
                        1e-8
                        //---rho_bar1---
                    ),
                    eps_constr(
                        //---eps_constr0---
                        1e-8
                        //---eps_constr1---
                    ),
                    xi_qn(
                        //---xi_qn0---
                        1e-4
                        //---xi_qn1---
                    ),
                    xi_pg(
                        //---xi_pg0---
                        1e-4
                        //---xi_pg1---
                    ),
                    xi_proj(
                        //---xi_proj0---
                        1e-4
                        //---xi_proj1---
                    ),
                    xi_tang(
                        //---xi_tang0---
                        1e-4
                        //---xi_tang1---
                    ),
                    xi_lmh(
                        //---xi_lmh0---
                        1e-4
                        //---xi_lmh1---
                    ),

                        //---xi_all0---
                        // None
                        //---xi_all1---

                    xi_lmg(
                        //---xi_lmg0---
                        1e4
                        //---xi_lmg1---
                    ),
                    xi_4(
                        //---xi_40---
                        2.
                        //---xi_41---
                    ),
                    rpred(
                        //---rpred0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---rpred1---
                    ),
                    PSchur_left_type(
                        //---PSchur_left_type0---
                        Operators::Identity
                        //---PSchur_left_type1---
                    ),
                    PSchur_right_type(
                        //---PSchur_right_type0---
                        Operators::Identity
                        //---PSchur_right_type1---
                    ),
                    augsys_iter_max(
                        //---augsys_iter_max0---
                        100
                        //---augsys_iter_max1---
                    ),
                    augsys_rst_freq(
                        //---augsys_rst_freq0---
                        0
                        //---augsys_rst_freq1---
                    ),
                    augsys_qn_iter(
                        //---augsys_qn_iter0---
                        0
                        //---augsys_qn_iter1---
                    ),
                    augsys_pg_iter(
                        //---augsys_pg_iter0---
                        0
                        //---augsys_pg_iter1---
                    ),
                    augsys_proj_iter(
                        //---augsys_proj_iter0---
                        0
                        //---augsys_proj_iter1---
                    ),
                    augsys_tang_iter(
                        //---augsys_tang_iter0---
                        0
                        //---augsys_tang_iter1---
                    ),
                    augsys_lmh_iter(
                        //---augsys_lmh_iter0---
                        0
                        //---augsys_lmh_iter1---
                    ),
                    augsys_qn_iter_total(
                        //---augsys_qn_iter_total0---
                        0
                        //---augsys_qn_iter_total1---
                    ),
                    augsys_pg_iter_total(
                        //---augsys_pg_iter_total0---
                        0
                        //---augsys_pg_iter_total1---
                    ),
                    augsys_proj_iter_total(
                        //---augsys_proj_iter_total0---
                        0
                        //---augsys_proj_iter_total1---
                    ),
                    augsys_tang_iter_total(
                        //---augsys_tang_iter_total0---
                        0
                        //---augsys_tang_iter_total1---
                    ),
                    augsys_lmh_iter_total(
                        //---augsys_lmh_iter_total0---
                        0
                        //---augsys_lmh_iter_total1---
                    ),
                    augsys_iter_total(
                        //---augsys_iter_total0---
                        0
                        //---augsys_iter_total1---
                    ),
                    augsys_qn_err(
                        //---augsys_qn_err0---
                        0.
                        //---augsys_qn_err1---
                    ),
                    augsys_pg_err(
                        //---augsys_pg_err0---
                        0.
                        //---augsys_pg_err1---
                    ),
                    augsys_proj_err(
                        //---augsys_proj_err0---
                        0.
                        //---augsys_proj_err1---
                    ),
                    augsys_tang_err(
                        //---augsys_tang_err0---
                        0.
                        //---augsys_tang_err1---
                    ),
                    augsys_lmh_err(
                        //---augsys_lmh_err0---
                        0.
                        //---augsys_lmh_err1---
                    ),
                    augsys_qn_err_target(
                        //---augsys_qn_err_target0---
                        0.
                        //---augsys_qn_err_target1---
                    ),
                    augsys_pg_err_target(
                        //---augsys_pg_err_target0---
                        0.
                        //---augsys_pg_err_target1---
                    ),
                    augsys_proj_err_target(
                        //---augsys_proj_err_target0---
                        0.
                        //---augsys_proj_err_target1---
                    ),
                    augsys_tang_err_target(
                        //---augsys_tang_err_target0---
                        0.
                        //---augsys_tang_err_target1---
                    ),
                    augsys_lmh_err_target(
                        //---augsys_lmh_err_target0---
                        0.
                        //---augsys_lmh_err_target1---
                    ),
                    g_x(
                        //---g_x0---
                        Y::init(y_user)
                        //---g_x1---
                    ),
                    norm_gxtyp(
                        //---norm_gxtyp0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---norm_gxtyp1---
                    ),
                    gpxdxn_p_gx(
                        //---gpxdxn_p_gx0---
                        Y::init(y_user)
                        //---gpxdxn_p_gx1---
                    ),
                    gpxdxt(
                        //---gpxdxt0---
                        Y::init(y_user)
                        //---gpxdxt1---
                    ),
                    norm_gpxdxnpgx(
                        //---norm_gpxdxnpgx0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---norm_gpxdxnpgx1---
                    ),
                    dx_n(
                        //---dx_n0---
                        X::init(x_user)
                        //---dx_n1---
                    ),
                    dx_ncp(
                        //---dx_ncp0---
                        X::init(x_user)
                        //---dx_ncp1---
                    ),
                    dx_t(
                        //---dx_t0---
                        X::init(x_user)
                        //---dx_t1---
                    ),
                    dx_t_uncorrected(
                        //---dx_t_uncorrected0---
                        X::init(x_user)
                        //---dx_t_uncorrected1---
                    ),
                    dx_tcp_uncorrected(
                        //---dx_tcp_uncorrected0---
                        X::init(x_user)
                        //---dx_tcp_uncorrected1---
                    ),
                    H_dxn(
                        //---H_dxn0---
                        X::init(x_user)
                        //---H_dxn1---
                    ),
                    W_gradpHdxn(
                        //---W_gradpHdxn0---
                        X::init(x_user)
                        //---W_gradpHdxn1---
                    ),
                    H_dxtuncorrected(
                        //---H_dxtuncorrected0---
                        X::init(x_user)
                        //---H_dxtuncorrected1---
                    ),
                    g_diag(
                        //---g_diag0---
                        FunctionDiagnostics::NoDiagnostics
                        //---g_diag1---
                    ),
                    y_diag(
                        //---y_diag0---
                        VectorSpaceDiagnostics::NoDiagnostics
                        //---y_diag1---
                    )
                {
                        Y::copy(y_user,y);
                }
            };
            
            // Check that we have a valid set of parameters.  
            static void check_(Messaging const & msg,t const & state) {
                   
                // Use this to build an error message
                std::stringstream ss;
                    
                    //---y_valid0---
                    // Any
                    //---y_valid1---
                    
                    //---dy_valid0---
                    // Any
                    //---dy_valid1---
                
                // Check that the fraction of the trust-region used for the
                // quasi-Newton step is between 0 and 1
                // is positive
                if(!(
                    //---zeta_valid0---
                    state.zeta > Real(0.) && state.zeta < Real(1.)
                    //---zeta_valid1---
                )) 
                    ss << "The fraction of the trust-region used for the "
                        "quasi-Newton step must lie in the interval (0,1): "
                        "zeta = " << state.zeta;
                
                // Check that the trust-region parameter that bounds the
                // error in the preduction reduction lies between 0 and 1-eta1.
                else if(!(
                    //---eta0_valid0---
                    state.eta0 > Real(0.) && state.eta0 < Real(1.)-state.eta1
                    //---eta0_valid1---
                )) 
                    ss << "The trust-region parameter that bounds the error "
                        "in the predicted reduction must lie in the interval "
                        "(0,1-eta1): eta0 = " << state.eta0;

                // Check that the augmented Lagrangian penalty parameter is
                // greater than or equal to 1
                else if(!(
                    //---rho_valid0---
                    state.rho >= Real(1.)
                    //---rho_valid1---
                ))
                    ss << "The augmented Lagrangian penalty parameter must be "
                        "greater than or equal to 1: rho = " << state.rho;

                // Check that the last penalty parameter is greater than or
                // equal to 1
                else if(!(
                    //---rho_old_valid0---
                    state.rho_old >= Real(1.)
                    //---rho_old_valid1---
                ))
                    ss << "The previous augmented Lagrangian penalty parameter"
                        "must be greater than or equal to 1: rho_old = "
                        << state.rho_old;

                // Check that the fixed increased to the augmented Lagrangian
                // penalty parameter is positive 
                else if(!(
                    //---rho_bar_valid0---
                    state.rho_bar > Real(0.)
                    //---rho_bar_valid1---
                ))
                    ss << "The fixed increase to the augmented Lagrangian "
                        "penalty paramter must be positive: rho_bar = " 
                        << state.rho_bar;

                // Check that the stopping tolerance for the norm of the
                // constraints is positive
                else if(!(
                    //---eps_constr_valid0---
                    state.eps_constr > Real(0.)
                    //---eps_constr_valid1---
                ))
                    ss << "The tolerance used in the norm of the constraints "
                        "stopping condition must be positive: eps_constr = "
                        << state.eps_constr;

                // Check that the quasi-Newton step inexactness tolerance lies 
                // in the interval (0,1) 
                else if(!(
                    //---xi_qn_valid0---
                    state.xi_qn > Real(0.) && state.xi_qn < Real(1.)
                    //---xi_qn_valid1---
                ))
                    ss << "The quasi-Newton step inexactness tolerance must "
                        "lie in the interval (0,1): xi_qn = " << state.xi_qn;
                
                // Check that the projected gradient inexactness tolerance lies 
                // in the interval (0,1) 
                else if(!(
                    //---xi_pg_valid0---
                    state.xi_pg > Real(0.) && state.xi_pg < Real(1.)
                    //---xi_pg_valid1---
                ))
                    ss << "The projected gradient inexactness tolerance must "
                        "lie in the interval (0,1): xi_pg = " << state.xi_pg;
                
                // Check that the nullspace projection inexactness tolerance
                // lies in the interval (0,1) 
                else if(!(
                    //---xi_proj_valid0---
                    state.xi_proj > Real(0.) && state.xi_proj < Real(1.)
                    //---xi_proj_valid1---
                ))
                    ss << "The nullspace projection inexactness tolerance must "
                        "lie in the interval (0,1): xi_proj = "<< state.xi_proj;
                
                // Check that the tangential step inexactness tolerance
                // lies in the interval (0,1) 
                else if(!(
                    //---xi_tang_valid0---
                    state.xi_tang > Real(0.) && state.xi_tang < Real(1.)
                    //---xi_tang_valid1---
                ))
                    ss << "The tangential step inexactness tolerance must "
                        "lie in the interval (0,1): xi_tang = "<< state.xi_tang;
                
                // Check that the equality multiplier inexactness tolerance
                // lies in the interval (0,1) 
                else if(!(
                    //---xi_lmh_valid0---
                    state.xi_lmh > Real(0.) && state.xi_lmh < Real(1.)
                    //---xi_lmh_valid1---
                ))
                    ss << "The equality multiplier inexactness tolerance must "
                        "lie in the interval (0,1): xi_lmh = " << state.xi_lmh;
                    
                    //---xi_all_valid0---
                    // state.xi_all > Real(0.) && state.xi_all < Real(1.)
                    //---xi_all_valid1---

                // Check that the absolute tolerance on the residual of the
                // equality multiplier solve is positive
                else if(!(
                    //---xi_lmg_valid0---
                    state.xi_lmg > Real(0.)
                    //---xi_lmg_valid1---
                ))
                    ss << "The equality multiplier residual tolerance must "
                        "be positive: xi_lmg = " << state.xi_lmg;

                // Check that the tolerance for the error acceptable in
                // the tangential step is greater than 1.
                else if(!(
                    //---xi_4_valid0---
                    state.xi_4 > Real(1.)
                    //---xi_4_valid1---
                ))
                    ss << "The tolerance on the acceptable error in the "
                        "tangential step must be greater than or equal to 1: "
                        "xi_4 = " << state.xi_4;
                    
                    //---rpred_valid0---
                    // Any
                    //---rpred_valid1---
                
                // Check that the left preconditioner for the augmented system
                // is either defined by the user or the identity.
                else if(!(
                    //---PSchur_left_type_valid0---
                    state.PSchur_left_type == Operators::Identity || 
                    state.PSchur_left_type == Operators::UserDefined
                    //---PSchur_left_type_valid1---
                ))
                    ss << "The left preconditioner for the augmented system "
                        "must be either user defined or the identity: "
                        "PSchur_left_type = "
                        << Operators::to_string(state.PSchur_left_type);
                
                // Check that the right preconditioner for the augmented system
                // is either defined by the user or the identity.
                else if(!(
                    //---PSchur_right_type_valid0---
                    state.PSchur_right_type == Operators::Identity || 
                    state.PSchur_right_type == Operators::UserDefined
                    //---PSchur_right_type_valid1---
                ))
                    ss << "The right preconditioner for the augmented system "
                        "must be either user defined or the identity: "
                        "PSchur_right_type = "
                        << Operators::to_string(state.PSchur_right_type);

                // Check that the number of iterations used when solving the
                // augmented system is positive
                else if(!(
                    //---augsys_iter_max_valid0---
                    state.augsys_iter_max > 0
                    //---augsys_iter_max_valid1---
                ))
                    ss << "The number of iterations used when solving the "
                        "augmented system must be positive: augsys_iter_max = "
                        << state.augsys_iter_max;
                    
                    //---augsys_rst_freq_valid0---
                    // Any
                    //---augsys_rst_freq_valid1---
                    
                    //---augsys_qn_iter_valid0---
                    // Any
                    //---augsys_qn_iter_valid1---
                    
                    //---augsys_pg_iter_valid0---
                    // Any
                    //---augsys_pg_iter_valid1---
                    
                    //---augsys_proj_iter_valid0---
                    // Any
                    //---augsys_proj_iter_valid1---
                    
                    //---augsys_tang_iter_valid0---
                    // Any
                    //---augsys_tang_iter_valid1---
                    
                    //---augsys_lmh_iter_valid0---
                    // Any
                    //---augsys_lmh_iter_valid1---
                    
                    //---augsys_qn_iter_total_valid0---
                    // Any
                    //---augsys_qn_iter_total_valid1---
                    
                    //---augsys_pg_iter_total_valid0---
                    // Any
                    //---augsys_pg_iter_total_valid1---
                    
                    //---augsys_proj_iter_total_valid0---
                    // Any
                    //---augsys_proj_iter_total_valid1---
                    
                    //---augsys_tang_iter_total_valid0---
                    // Any
                    //---augsys_tang_iter_total_valid1---
                    
                    //---augsys_lmh_iter_total_valid0---
                    // Any
                    //---augsys_lmh_iter_total_valid1---
                    
                    //---augsys_iter_total_valid0---
                    // Any
                    //---augsys_iter_total_valid1---
                    
                    //---augsys_qn_err_valid0---
                    // Any
                    //---augsys_qn_err_valid1---
                    
                    //---augsys_pg_err_valid0---
                    // Any
                    //---augsys_pg_err_valid1---
                    
                    //---augsys_proj_err_valid0---
                    // Any
                    //---augsys_proj_err_valid1---
                    
                    //---augsys_tang_err_valid0---
                    // Any
                    //---augsys_tang_err_valid1---
                    
                    //---augsys_lmh_err_valid0---
                    // Any
                    //---augsys_lmh_err_valid1---
                    
                    //---augsys_qn_err_target_valid0---
                    // Any
                    //---augsys_qn_err_target_valid1---
                    
                    //---augsys_pg_err_target_valid0---
                    // Any
                    //---augsys_pg_err_target_valid1---
                    
                    //---augsys_proj_err_target_valid0---
                    // Any
                    //---augsys_proj_err_target_valid1---
                    
                    //---augsys_tang_err_target_valid0---
                    // Any
                    //---augsys_tang_err_target_valid1---
                    
                    //---augsys_lmh_err_target_valid0---
                    // Any
                    //---augsys_lmh_err_target_valid1---
                    
                    //---g_x_valid0---
                    // Any
                    //---g_x_valid1---

                // Check that the norm of a typical constraint is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(!(
                    //---norm_gxtyp_valid0---
                    state.norm_gxtyp >= Real(0.)
                    || (state.iter==1 && state.norm_gxtyp!=state.norm_gxtyp)
                    //---norm_gxtyp_valid1---
                )) 
                    ss << "The norm of a typical constraint must be "
                        "nonnegative: norm_gxtyp = " << state.norm_gxtyp; 
                    
                    //---gpxdxn_p_gx_valid0---
                    // Any
                    //---gpxdxn_p_gx_valid1---
                    
                    //---gpxdxt_valid0---
                    // Any
                    //---gpxdxt_valid1---
                    
                    //---norm_gpxdxnpgx_valid0---
                    // Any
                    //---norm_gpxdxnpgx_valid1---
                    
                    //---dx_n_valid0---
                    // Any
                    //---dx_n_valid1---
                    
                    //---dx_ncp_valid0---
                    // Any
                    //---dx_ncp_valid1---
                    
                    //---dx_t_valid0---
                    // Any
                    //---dx_t_valid1---
                    
                    //---dx_t_uncorrected_valid0---
                    // Any
                    //---dx_t_uncorrected_valid1---
                    
                    //---dx_tcp_uncorrected_valid0---
                    // Any
                    //---dx_tcp_uncorrected_valid1---
                    
                    //---H_dxn_valid0---
                    // Any
                    //---H_dxn_valid1---
                    
                    //---W_gradpHdxn_valid0---
                    // Any
                    //---W_gradpHdxn_valid1---
                    
                    //---H_dxtuncorrected_valid0---
                    // Any
                    //---H_dxtuncorrected_valid1---
                    
                    //---g_diag_valid0---
                    // Any
                    //---g_diag_valid1---
                    
                    //---y_diag_valid0---
                    // Any
                    //---y_diag_valid1---

                // If there's an error, print it
                if(ss.str()!="") msg.error(ss.str());
            }
            static void check(Messaging const & msg,t const & state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {
            // Disallow constructors
            NO_CONSTRUCTORS(Restart)
       
            // Create some type shortcuts 
            typedef typename RestartPackage <Real>::t Reals;
            typedef typename RestartPackage <Natural>::t Naturals;
            typedef typename RestartPackage <std::string>::t Params;
            typedef typename RestartPackage <X_Vector>::t X_Vectors;
            typedef typename RestartPackage <Y_Vector>::t Y_Vectors;

            // Checks whether we have a valid real 
            static bool is_real(
                typename RestartPackage <Real>::tuple const & item
            ){
                if( Unconstrained <Real,XX>::Restart::is_real(item) ||
                    item.first == "zeta" ||
                    item.first == "eta0" ||
                    item.first == "rho" ||
                    item.first == "rho_old" ||
                    item.first == "rho_bar" ||
                    item.first == "eps_constr" ||
                    item.first == "xi_qn" || 
                    item.first == "xi_pg" ||
                    item.first == "xi_proj" ||
                    item.first == "xi_tang" ||
                    item.first == "xi_lmh" ||
                    item.first == "xi_lmg" ||
                    item.first == "xi_4" ||
                    item.first == "rpred" ||
                    item.first == "norm_gxtyp" ||
                    item.first == "norm_gpxdxnpgx" ||
                    item.first == "augsys_qn_err" ||
                    item.first == "augsys_pg_err" ||
                    item.first == "augsys_proj_err" ||
                    item.first == "augsys_tang_err" ||
                    item.first == "augsys_lmh_err" ||
                    item.first == "augsys_qn_err_target" ||
                    item.first == "augsys_pg_err_target" ||
                    item.first == "augsys_proj_err_target" ||
                    item.first == "augsys_tang_err_target" ||
                    item.first == "augsys_lmh_err_target"
                )
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid natural number
            static bool is_nat(
                typename RestartPackage <Natural>::tuple const & item
            ) {
                if( Unconstrained <Real,XX>::Restart::is_nat(item) ||
                    item.first == "augsys_iter_max" ||
                    item.first == "augsys_rst_freq" ||
                    item.first == "augsys_qn_iter" ||
                    item.first == "augsys_pg_iter" ||
                    item.first == "augsys_proj_iter" ||
                    item.first == "augsys_tang_iter" ||
                    item.first == "augsys_lmh_iter" ||
                    item.first == "augsys_qn_iter_total" ||
                    item.first == "augsys_pg_iter_total" ||
                    item.first == "augsys_proj_iter_total" ||
                    item.first == "augsys_tang_iter_total" ||
                    item.first == "augsys_lmh_iter_total" ||
                    item.first == "augsys_iter_total"
                )
                    return true;
                else
                    return false;
            }
           
            // Checks whether we have a valid parameter 
            static bool is_param(
                typename RestartPackage <std::string>::tuple const & item
            ){
                if( Unconstrained <Real,XX>::Restart::is_param(item) ||
                    (item.first=="PSchur_left_type" &&
                        Operators::is_valid(item.second)) ||
                    (item.first=="PSchur_right_type" &&
                        Operators::is_valid(item.second)) ||
                    (item.first=="g_diag" &&
                        FunctionDiagnostics::is_valid(item.second)) ||
                    (item.first=="y_diag" &&
                        VectorSpaceDiagnostics::is_valid(item.second))
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid variable
            static bool is_x (
                typename RestartPackage <X_Vector>::tuple const & item
            ) {
                if( Unconstrained <Real,XX>::Restart::is_x(item) ||
                    item.first == "dx_n" ||
                    item.first == "dx_ncp" ||
                    item.first == "dx_t" ||
                    item.first == "dx_t_uncorrected" ||
                    item.first == "dx_tcp_uncorrected" ||
                    item.first == "H_dxn" ||
                    item.first == "W_gradpHdxn" ||
                    item.first == "H_dxtuncorrected" 
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid equality multiplier 
            static bool is_y(
                typename RestartPackage <Y_Vector>::tuple const & item
            ) {
                if( item.first == "y" ||
                    item.first == "dy" ||
                    item.first == "g_x" ||
                    item.first == "gpxdxn_p_gx" ||
                    item.first == "gpxdxt"
                ) 
                    return true;
                else
                    return false;
            }

            // Checks whether we have valid labels
            static void checkItems(
                Messaging const & msg,
                Reals const & reals,
                Naturals const & nats,
                Params const & params,
                X_Vectors const & xs,
                Y_Vectors const & ys
            ) {
                Utility::checkItems <Real> (
                    msg,is_real,reals," real name: ");
                Utility::checkItems <Natural> (
                    msg,is_nat,nats," natural name: ");
                Utility::checkItems <std::string> (
                    msg,is_param,params," paramater: ");
                Utility::checkItems <X_Vector> (
                    msg,is_x,xs," variable name: ");
                Utility::checkItems <Y_Vector> (
                    msg,is_y,ys," equality multiplier name: ");
            }
            
            // Copy out all equality multipliers 
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs,
                Y_Vectors & ys
            ) {
                ys.emplace_back("y",std::move(state.y));
                ys.emplace_back("dy",std::move(state.dy));
                ys.emplace_back("g_x",std::move(state.g_x));
                ys.emplace_back("gpxdxn_p_gx",std::move(state.gpxdxn_p_gx));
                ys.emplace_back("gpxdxt",std::move(state.gpxdxt));
                
                xs.emplace_back("dx_n",std::move(state.dx_n));
                xs.emplace_back("dx_ncp",std::move(state.dx_ncp));
                xs.emplace_back("dx_t",std::move(state.dx_t));
                xs.emplace_back("dx_t_uncorrected",
                    std::move(state.dx_t_uncorrected));
                xs.emplace_back("dx_tcp_uncorrected",
                    std::move(state.dx_tcp_uncorrected));
                xs.emplace_back("H_dxn",std::move(state.H_dxn));
                xs.emplace_back("W_gradpHdxn",std::move(state.W_gradpHdxn));
                xs.emplace_back("H_dxtuncorrected",
                    std::move(state.H_dxtuncorrected));
            }

            // Copy out all the scalar information
            static void stateToScalars(
                typename State::t & state,
                Reals & reals,
                Naturals & nats,
                Params & params
            ) { 
                // Copy in all the real numbers
                reals.emplace_back("zeta",std::move(state.zeta));
                reals.emplace_back("eta0",std::move(state.eta0));
                reals.emplace_back("rho",std::move(state.rho));
                reals.emplace_back("rho_old",std::move(state.rho_old));
                reals.emplace_back("rho_bar",std::move(state.rho_bar));
                reals.emplace_back("eps_constr",std::move(state.eps_constr));
                reals.emplace_back("xi_qn",std::move(state.xi_qn));
                reals.emplace_back("xi_pg",std::move(state.xi_pg));
                reals.emplace_back("xi_proj",std::move(state.xi_proj));
                reals.emplace_back("xi_tang",std::move(state.xi_tang));
                reals.emplace_back("xi_lmh",std::move(state.xi_lmh));
                reals.emplace_back("xi_lmg",std::move(state.xi_lmg));
                reals.emplace_back("xi_4",std::move(state.xi_4));
                reals.emplace_back("rpred",std::move(state.rpred));
                reals.emplace_back("norm_gxtyp",std::move(state.norm_gxtyp));
                reals.emplace_back("norm_gpxdxnpgx",
                    std::move(state.norm_gpxdxnpgx));
                reals.emplace_back("augsys_qn_err",
                    std::move(state.augsys_qn_err));
                reals.emplace_back("augsys_pg_err",
                    std::move(state.augsys_pg_err));
                reals.emplace_back("augsys_proj_err",
                    std::move(state.augsys_proj_err));
                reals.emplace_back("augsys_tang_err",
                    std::move(state.augsys_tang_err));
                reals.emplace_back("augsys_lmh_err",
                    std::move(state.augsys_lmh_err));
                reals.emplace_back("augsys_qn_err_target",
                    std::move(state.augsys_qn_err_target));
                reals.emplace_back("augsys_pg_err_target",
                    std::move(state.augsys_pg_err_target));
                reals.emplace_back("augsys_proj_err_target",
                    std::move(state.augsys_proj_err_target));
                reals.emplace_back("augsys_tang_err_target",
                    std::move(state.augsys_tang_err_target));
                reals.emplace_back("augsys_lmh_err_target",
                    std::move(state.augsys_lmh_err_target));

                // Copy in all the natural numbers
                nats.emplace_back("augsys_iter_max",
                    std::move(state.augsys_iter_max));
                nats.emplace_back("augsys_rst_freq",
                    std::move(state.augsys_rst_freq));
                nats.emplace_back("augsys_qn_iter",
                    std::move(state.augsys_qn_iter));
                nats.emplace_back("augsys_pg_iter",
                    std::move(state.augsys_pg_iter));
                nats.emplace_back("augsys_proj_iter",
                    std::move(state.augsys_proj_iter));
                nats.emplace_back("augsys_tang_iter",
                    std::move(state.augsys_tang_iter));
                nats.emplace_back("augsys_lmh_iter",
                    std::move(state.augsys_lmh_iter));
                nats.emplace_back("augsys_qn_iter_total",
                    std::move(state.augsys_qn_iter_total));
                nats.emplace_back("augsys_pg_iter_total",
                    std::move(state.augsys_pg_iter_total));
                nats.emplace_back("augsys_proj_iter_total",
                    std::move(state.augsys_proj_iter_total));
                nats.emplace_back("augsys_tang_iter_total",
                    std::move(state.augsys_tang_iter_total));
                nats.emplace_back("augsys_lmh_iter_total",
                    std::move(state.augsys_lmh_iter_total));
                nats.emplace_back("augsys_iter_total",
                    std::move(state.augsys_iter_total));

                // Copy in all the parameters
                params.emplace_back("PSchur_left_type",
                    Operators::to_string(state.PSchur_left_type));
                params.emplace_back("PSchur_right_type",
                    Operators::to_string(state.PSchur_right_type));
                params.emplace_back("g_diag",
                    FunctionDiagnostics::to_string(state.g_diag));
                params.emplace_back("y_diag",
                    VectorSpaceDiagnostics::to_string(state.y_diag));
            }
            
            // Copy in all equality multipliers 
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys
            ) {
                for(typename Y_Vectors::iterator item = ys.begin();
                    item!=ys.end();
                    item++
                ){
                    if(item->first=="y")
                        state.y = std::move(item->second);
                    else if(item->first=="dy")
                        state.dy = std::move(item->second);
                    else if(item->first=="g_x")
                        state.g_x = std::move(item->second);
                    else if(item->first=="gpxdxn_p_gx")
                        state.gpxdxn_p_gx = std::move(item->second);
                    else if(item->first=="gpxdxt")
                        state.gpxdxt = std::move(item->second);
                }

                for(typename X_Vectors::iterator item = xs.begin();
                    item!=xs.end();
                    item++
                ){
                    if(item->first=="dx_n")
                        state.dx_n = std::move(item->second);
                    else if(item->first=="dx_ncp")
                        state.dx_ncp = std::move(item->second);
                    else if(item->first=="dx_t")
                        state.dx_t = std::move(item->second);
                    else if(item->first=="dx_t_uncorrected")
                        state.dx_t_uncorrected = std::move(item->second);
                    else if(item->first=="dx_tcp_uncorrected")
                        state.dx_tcp_uncorrected = std::move(item->second);
                    else if(item->first=="H_dxn")
                        state.H_dxn = std::move(item->second);
                    else if(item->first=="W_gradpHdxn")
                        state.W_gradpHdxn = std::move(item->second);
                    else if(item->first=="H_dxtuncorrected")
                        state.H_dxtuncorrected = std::move(item->second);
                }
            }
            
            // Copy in all the scalar information
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Naturals & nats,
                Params & params
            ) { 
                // Copy in any reals 
                for(typename Reals::iterator item = reals.begin();
                    item!=reals.end();
                    item++
                ){
                    if(item->first=="zeta")
                        state.zeta=std::move(item->second);
                    else if(item->first=="eta0")
                        state.eta0=std::move(item->second);
                    else if(item->first=="rho")
                        state.rho=std::move(item->second);
                    else if(item->first=="rho_old")
                        state.rho_old=std::move(item->second);
                    else if(item->first=="rho_bar")
                        state.rho_bar=std::move(item->second);
                    else if(item->first=="eps_constr")
                        state.eps_constr=std::move(item->second);
                    else if(item->first=="xi_qn")
                        state.xi_qn=std::move(item->second);
                    else if(item->first=="xi_pg")
                        state.xi_pg=std::move(item->second);
                    else if(item->first=="xi_proj")
                        state.xi_proj=std::move(item->second);
                    else if(item->first=="xi_tang")
                        state.xi_tang=std::move(item->second);
                    else if(item->first=="xi_lmh")
                        state.xi_lmh=std::move(item->second);
                    else if(item->first=="xi_lmg")
                        state.xi_lmg=std::move(item->second);
                    else if(item->first=="xi_4")
                        state.xi_4=std::move(item->second);
                    else if(item->first=="rpred")
                        state.rpred=std::move(item->second);
                    else if(item->first=="norm_gxtyp")
                        state.norm_gxtyp=std::move(item->second);
                    else if(item->first=="norm_gpxdxnpgx")
                        state.norm_gpxdxnpgx=std::move(item->second);
                    else if(item->first=="augsys_qn_err")
                        state.augsys_qn_err=std::move(item->second);
                    else if(item->first=="augsys_pg_err")
                        state.augsys_pg_err=std::move(item->second);
                    else if(item->first=="augsys_proj_err")
                        state.augsys_proj_err=std::move(item->second);
                    else if(item->first=="augsys_tang_err")
                        state.augsys_tang_err=std::move(item->second);
                    else if(item->first=="augsys_lmh_err")
                        state.augsys_lmh_err=std::move(item->second);
                    else if(item->first=="augsys_qn_err_target")
                        state.augsys_qn_err_target=std::move(item->second);
                    else if(item->first=="augsys_pg_err_target")
                        state.augsys_pg_err_target=std::move(item->second);
                    else if(item->first=="augsys_proj_err_target")
                        state.augsys_proj_err_target=std::move(item->second);
                    else if(item->first=="augsys_tang_err_target")
                        state.augsys_tang_err_target=std::move(item->second);
                    else if(item->first=="augsys_lmh_err_target")
                        state.augsys_lmh_err_target=std::move(item->second);
                }
                
                // Next, copy in any naturals
                for(typename Naturals::iterator item = nats.begin();
                    item!=nats.end();
                    item++
                ){
                    if(item->first=="augsys_iter_max")
                        state.augsys_iter_max=std::move(item->second);
                    else if(item->first=="augsys_rst_freq")
                        state.augsys_rst_freq=std::move(item->second);
                    else if(item->first=="augsys_qn_iter")
                        state.augsys_qn_iter=std::move(item->second);
                    else if(item->first=="augsys_pg_iter")
                        state.augsys_pg_iter=std::move(item->second);
                    else if(item->first=="augsys_proj_iter")
                        state.augsys_proj_iter=std::move(item->second);
                    else if(item->first=="augsys_tang_iter")
                        state.augsys_tang_iter=std::move(item->second);
                    else if(item->first=="augsys_lmh_iter")
                        state.augsys_lmh_iter=std::move(item->second);
                    else if(item->first=="augsys_qn_iter_total")
                        state.augsys_qn_iter_total=std::move(item->second);
                    else if(item->first=="augsys_pg_iter_total")
                        state.augsys_pg_iter_total=std::move(item->second);
                    else if(item->first=="augsys_proj_iter_total")
                        state.augsys_proj_iter_total=std::move(item->second);
                    else if(item->first=="augsys_tang_iter_total")
                        state.augsys_tang_iter_total=std::move(item->second);
                    else if(item->first=="augsys_lmh_iter_total")
                        state.augsys_lmh_iter_total=std::move(item->second);
                    else if(item->first=="augsys_iter_total")
                        state.augsys_iter_total=std::move(item->second);
                }
                
                // Next, copy in any parameters 
                for(typename Params::iterator item = params.begin();
                    item!=params.end();
                    item++
                ){
                    if(item->first=="PSchur_left_type")
                        state.PSchur_left_type
                            = Operators::from_string(item->second);
                    else if(item->first=="PSchur_right_type")
                        state.PSchur_right_type
                            =Operators::from_string(item->second);
                    else if(item->first=="g_diag")
                        state.g_diag=FunctionDiagnostics::from_string(
                            item->second);
                    else if(item->first=="y_diag")
                        state.y_diag=VectorSpaceDiagnostics::from_string(
                            item->second);
                }
            }

            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Reals & reals,
                Naturals & nats,
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
                Naturals & nats,
                Params & params
            ) {

                // Check the user input 
                checkItems(msg,reals,nats,params,xs,ys);

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
            NO_CONSTRUCTORS(Functions)
            
            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t)

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
                NO_COPY_ASSIGNMENT(EqualityModifications)

            private:
                // Underlying modification.  This takes control of the memory
                std::unique_ptr <
                    Optizelle::ScalarValuedFunctionModifications <Real,XX>
                > f_mod;

                // Equality constraint.
                Optizelle::VectorValuedFunction <Real,XX,YY> const & g;

                // Reference to equality multiplier
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
                    typename Functions::t const & fns,
                    typename State::t const & state,
                    std::unique_ptr <
                        ScalarValuedFunctionModifications<Real,XX>
                    > && f_mod_
                ) : f_mod(std::move(f_mod_)),
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
                        g.eval(x,g_x);
                    
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
                void eval(Y_Vector const & dy,Y_Vector & result) const{
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

                // Check that all functions are defined 
                check(msg,fns);
                
                // Modify the objective 
                fns.f_mod.reset(new EqualityModifications(
                    fns,state,std::move(fns.f_mod)));
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
        struct Diagnostics {
            // Disallow constructors
            NO_CONSTRUCTORS(Diagnostics)

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { 
                // Create some shortcuts
                Natural const & msg_level = state.msg_level; 

                // Norm of the constrained 
                out.emplace_back(Utility::atos("||g(x)||"));
                
                // More detailed information
                if(msg_level>=2) {
                    // Trust-region information
                    out.emplace_back(Utility::atos("ared"));
                    out.emplace_back(Utility::atos("pred"));
                    out.emplace_back(Utility::atos("ared/pred"));
                    out.emplace_back(Utility::atos("delta"));
                       
                    // Krylov method information
                    out.emplace_back(Utility::atos("kry_iter"));
                    out.emplace_back(Utility::atos("kry_err"));
                    out.emplace_back(Utility::atos("kry_why"));
                }

                // Even more detail
                if(msg_level>=3) {
                    // Size of the normal and tangential steps
                    out.emplace_back(Utility::atos("|| dx_n ||"));
                    out.emplace_back(Utility::atos("|| dx_t ||"));

                    // Total number of Krylov iterations
                    out.emplace_back(Utility::atos("kry_itr_tot"));

                    // Augmented system solves 
                    out.emplace_back(Utility::atos("qn_iter"));
                    out.emplace_back(Utility::atos("qn_iter_tot"));
                    out.emplace_back(Utility::atos("qn_err"));
                    out.emplace_back(Utility::atos("qn_err_trg"));
                    
                    out.emplace_back(Utility::atos("pg_iter"));
                    out.emplace_back(Utility::atos("pg_iter_tot"));
                    out.emplace_back(Utility::atos("pg_err"));
                    out.emplace_back(Utility::atos("pg_err_trg"));
                    
                    out.emplace_back(Utility::atos("pr_iter"));
                    out.emplace_back(Utility::atos("pr_iter_tot"));
                    out.emplace_back(Utility::atos("pr_err"));
                    out.emplace_back(Utility::atos("pr_err_trg"));
                    
                    out.emplace_back(Utility::atos("tg_iter"));
                    out.emplace_back(Utility::atos("tg_iter_tot"));
                    out.emplace_back(Utility::atos("tg_err"));
                    out.emplace_back(Utility::atos("tg_err_trg"));
                    
                    out.emplace_back(Utility::atos("lm_iter"));
                    out.emplace_back(Utility::atos("lm_iter_tot"));
                    out.emplace_back(Utility::atos("lm_err"));
                    out.emplace_back(Utility::atos("lm_err_trg"));
                    
                    out.emplace_back(Utility::atos("aug_itr_tot"));
                }
            }
            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getStateHeader_(
                    state,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getStateHeader_(
                    state,out);
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
                Natural const & krylov_iter_total=state.krylov_iter_total;
                Real const & krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t const & krylov_stop=state.krylov_stop;
                Real const & pred = state.pred;
                Real const & ared = state.ared;
                Real const & delta = state.delta;
                Natural const & msg_level = state.msg_level;
                auto const & algorithm_class = state.algorithm_class;

                auto const & dx_n = state.dx_n;
                auto const & dx_t = state.dx_t;

                auto const & augsys_qn_iter = state.augsys_qn_iter;
                auto const & augsys_qn_iter_total=state.augsys_qn_iter_total;
                auto const & augsys_qn_err = state.augsys_qn_err;
                auto const & augsys_qn_err_target= state.augsys_qn_err_target;

                auto const & augsys_pg_iter = state.augsys_pg_iter;
                auto const & augsys_pg_iter_total=state.augsys_pg_iter_total;
                auto const & augsys_pg_err = state.augsys_pg_err;
                auto const & augsys_pg_err_target= state.augsys_pg_err_target;

                auto const & augsys_proj_iter = state.augsys_proj_iter;
                auto const & augsys_proj_iter_total =
                    state.augsys_proj_iter_total;
                auto const & augsys_proj_err = state.augsys_proj_err;
                auto const & augsys_proj_err_target =
                    state.augsys_proj_err_target;

                auto const & augsys_tang_iter = state.augsys_tang_iter;
                auto const & augsys_tang_iter_total =
                    state.augsys_tang_iter_total;
                auto const & augsys_tang_err = state.augsys_tang_err;
                auto const & augsys_tang_err_target =
                    state.augsys_tang_err_target;

                auto const & augsys_lmh_iter = state.augsys_lmh_iter;
                auto const & augsys_lmh_iter_total =
                    state.augsys_lmh_iter_total;
                auto const & augsys_lmh_err = state.augsys_lmh_err;
                auto const & augsys_lmh_err_target =
                    state.augsys_lmh_err_target;
                
                auto const & augsys_iter_total = state.augsys_iter_total;

                // Figure out if we're at the absolute beginning of the
                // optimization.
                bool opt_begin = Utility::is_opt_begin <EqualityConstrained> (
                    state);

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Norm of the gradient 
                Real norm_gx = sqrt(Y::innr(g_x,g_x));
                out.emplace_back(Utility::atos(norm_gx));
                    
                // More detailed information
                if(msg_level >=2) {
                    // Actual vs. predicted reduction 
                    if(!opt_begin) {
                        out.emplace_back(Utility::atos(ared));
                        out.emplace_back(Utility::atos(pred));
                        out.emplace_back(Utility::atos(ared/pred));
                        out.emplace_back(Utility::atos(delta));
                    } else 
                        for(Natural i=0;i<4;i++)
                            out.emplace_back(Utility::blankSeparator);
                    
                    // Krylov method information
                    if(!opt_begin) {
                        out.emplace_back(Utility::atos(krylov_iter));
                        out.emplace_back(Utility::atos(krylov_rel_err));
                        out.emplace_back(Utility::atos(krylov_stop));
                    } else 
                        for(Natural i=0;i<3;i++)
                            out.emplace_back(Utility::blankSeparator);
                }
                
                // Even more detail
                if(msg_level >=3) {
                    if(!opt_begin) {
                        // Size of the normal and tangential steps
                        auto norm_dxn = std::sqrt(X::innr(dx_n,dx_n));
                        auto norm_dxt = std::sqrt(X::innr(dx_t,dx_t));
                        out.emplace_back(Utility::atos(norm_dxn));
                        out.emplace_back(Utility::atos(norm_dxt));

                        // Total number of Krylov iterations
                        out.emplace_back(Utility::atos(krylov_iter_total));

                        // Augmented system solves 
                        out.emplace_back(Utility::atos(augsys_qn_iter));
                        out.emplace_back(Utility::atos(augsys_qn_iter_total));
                        out.emplace_back(Utility::atos(augsys_qn_err));
                        out.emplace_back(Utility::atos(augsys_qn_err_target));
                        
                        out.emplace_back(Utility::atos(augsys_pg_iter));
                        out.emplace_back(Utility::atos(augsys_pg_iter_total));
                        out.emplace_back(Utility::atos(augsys_pg_err));
                        out.emplace_back(Utility::atos(augsys_pg_err_target));
                        
                        out.emplace_back(Utility::atos(augsys_proj_iter));
                        out.emplace_back(Utility::atos(augsys_proj_iter_total));
                        out.emplace_back(Utility::atos(augsys_proj_err));
                        out.emplace_back(
                            Utility::atos(augsys_proj_err_target));
                        
                        out.emplace_back(Utility::atos(augsys_tang_iter));
                        out.emplace_back(Utility::atos(augsys_tang_iter_total));
                        out.emplace_back(Utility::atos(augsys_tang_err));
                        out.emplace_back(
                            Utility::atos(augsys_tang_err_target));
                        
                        out.emplace_back(Utility::atos(augsys_lmh_iter));
                        out.emplace_back(Utility::atos(augsys_lmh_iter_total));
                        out.emplace_back(Utility::atos(augsys_lmh_err));
                        out.emplace_back(Utility::atos(augsys_lmh_err_target));
                        
                        out.emplace_back(Utility::atos(augsys_iter_total));
                    } else 
                        for(Natural i=0;i<24;i++)
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
                Unconstrained <Real,XX>::Diagnostics
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics
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
                Unconstrained <Real,XX>::Diagnostics::getKrylovHeader_(
                    state,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getKrylovHeader_(
                    state,out);
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
                Unconstrained <Real,XX>::Diagnostics::getKrylov_(
                    state,blank,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getKrylov_(
                    state,blank,out);
            }
           
            // Runs the specified function diagnostics 
            static void checkFunctions_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);
                X_Vector const & x=state.x;
                Y_Vector const & y=state.y;
                FunctionDiagnostics::t const & g_diag = state.g_diag;
                
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 
                Y_Vector dy(Y::init(y));
                    Y::rand(dy);

                // Run the diagnostics
                switch(g_diag) {
                    case FunctionDiagnostics::FirstOrder:
                        msg.print("Diagnostics on the function g");
                        Optizelle::Diagnostics::derivativeCheck(
                            msg,g,x,dx,dy,"g");
                        Optizelle::Diagnostics::derivativeAdjointCheck(
                            msg,g,x,dx,dy,"g");
                        msg.print("");
                        break;
                    case FunctionDiagnostics::SecondOrder:
                        msg.print("Diagnostics on the function g");
                        Optizelle::Diagnostics::derivativeCheck(
                            msg,g,x,dx,dy,"g");
                        Optizelle::Diagnostics::derivativeAdjointCheck(
                            msg,g,x,dx,dy,"g");
                        Optizelle::Diagnostics::secondDerivativeCheck(
                            msg,g,x,dx,dy,"g");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified function diagnostics 
            static void checkFunctions(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkFunctions_(
                    msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics::checkFunctions_(
                    msg,fns,state);
            }
            
            // Sets up the equality Hessian operator 
            struct EqualityHessianOperator : public Operator <Real,XX,XX> {
            private:
                // Store the equality constraints 
                VectorValuedFunction <Real,XX,YY> const & g;

                // Current iterate
                X_Vector const & x; 

                // Equality modifications 
                typename Functions::EqualityModifications g_mod;

            public:
                // Remove some constructors
                NO_COPY_ASSIGNMENT(EqualityHessianOperator);

                // Take in the objective and the base point during construction 
                EqualityHessianOperator(
                    typename Functions::t const & fns,
                    typename State::t const & state
                ) :
                    g(*(fns.g)),
                    x(state.x),
                    g_mod(
                        fns,
                        state,
                        std::unique_ptr <
                            ScalarValuedFunctionModifications <Real,XX>
                        > (new ScalarValuedFunctionModifications <Real,XX> ())
                    )
                {}

                // Basic application
                void eval(X_Vector const & dx,X_Vector & result)
                    const
                {
                    // Grab the zero vector in the X space
                    X_Vector zero(X::init(x));
                    X::zero(zero);

                    // Add the equality constraint's contribution to the
                    // Hessian-vector product
                    g_mod.hessvec_step(x,dx,zero,result);
                }
            };
           
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                X_Vector const & x=state.x;
                FunctionDiagnostics::t const & L_diag=state.L_diag;
                
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 
                X_Vector dxx(X::init(x));
                    X::rand(dxx); 

                // Create the equality Hessian operator
                EqualityHessianOperator L(fns,state);

                // Run the diagnostics
                switch(L_diag) {
                    case FunctionDiagnostics::SecondOrder:
                        msg.print("Diagnostics on the contribution of g to "
                            "the Lagrangian");
                        Optizelle::Diagnostics::operatorSymmetryCheck <Real,XX>(
                            msg,L,dx,dxx,"(g''(x).)*y");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkLagrangian_(
                    msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics::checkLagrangian_(
                    msg,fns,state);
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                VectorSpaceDiagnostics::t const & y_diag=state.y_diag;
                Y_Vector const & y=state.y;
               
                // Create some random directions for these tests
                Y_Vector dy(Y::init(y));
                    Y::rand(dy); 

                // Run the diagnostics
                switch(y_diag) {
                    case VectorSpaceDiagnostics::Basic:
                        msg.print("Diagnostics on the vector-space Y");
                        Optizelle::Diagnostics::zero_innr <Real,YY> (msg,y,"Y");
                        Optizelle::Diagnostics::copy_axpy_innr <Real,YY> (
                            msg,dy,"Y");
                        Optizelle::Diagnostics::copy_scal_innr <Real,YY> (
                            msg,dy,"Y");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
            }
        };

        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms)

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
                void eval(
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
                void eval(
                    const XxY_Vector & dx_dy,
                    XxY_Vector & result
                ) const{
                    // PH_x dx
                    PH_x.eval(dx_dy.first,result.first);
                    
                    // PH_y dy
                    PH_y.eval(dx_dy.second,result.second);
                }
            };

            // Sets the tolerances for the quasi-normal Newton solve
            struct QNManipulator : public GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public:
                // Disallow constructors
                NO_DEFAULT_COPY_ASSIGNMENT(QNManipulator)

                // Grab the states and fns on construction
                explicit QNManipulator(
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}

                // Evalulate the manipulator
                void eval(
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    auto const & absrel = *(fns.absrel);
                    Real const & xi_qn = state.xi_qn;
                    Real const & norm_gxtyp = state.norm_gxtyp;
                    Real const & eps_constr= state.eps_constr;
                    Real & augsys_qn_err_target = state.augsys_qn_err_target;

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
                    if(norm_gpxdxncp_p_g < eps_constr * absrel(norm_gxtyp))
                        eps=Real(1.);

                    // Save this desired error
                    augsys_qn_err_target=eps;
                }
            };

            // Finds the quasi-normal step
            static void quasinormalStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,YY> const & g=*(fns.g);
                auto const & safeguard = *(fns.safeguard); 
                auto const & absrel = *(fns.absrel);
                X_Vector const & x=state.x;
                Y_Vector const & g_x=state.g_x;
                Natural const & augsys_iter_max=state.augsys_iter_max;
                Natural const & augsys_rst_freq=state.augsys_rst_freq;
                Real const & delta = state.delta;
                Real const & zeta = state.zeta;
                Real const & norm_gxtyp = state.norm_gxtyp;
                Real const & eps_constr = state.eps_constr;
                X_Vector & dx_ncp=state.dx_ncp;
                X_Vector & dx_n=state.dx_n;
                Real & augsys_qn_err = state.augsys_qn_err;
                Natural & augsys_qn_iter = state.augsys_qn_iter;
                Natural & augsys_qn_iter_total = state.augsys_qn_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;
                auto & alpha_x_qn = state.alpha_x_qn;

                // If we're already feasible, don't even bother with the
                // quasi-Newton step.  In fact, if g(x)=0, the equation for
                // the Cauchy point divides by zero, which causes all sorts
                // of headaches later on.
                Real norm_gx = sqrt(X::innr(g_x,g_x));
                if(norm_gx < eps_constr * absrel(norm_gxtyp)) {
                    X::zero(dx_ncp);
                    X::zero(dx_n);
                    return;
                }

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

                // Safeguard the Cauchy point
                auto zero = X::init(x);
                X::zero(zero);
                alpha_x_qn = std::min(safeguard(zero,dx_ncp,zeta),Real(1.));

                // If || dx_ncp || >= zeta delta, scale it back to zeta
                // delta and return.  Alternatively, if the safeguard truncates
                // things, scale things back and return.
                Real norm_dxncp = sqrt(X::innr(dx_ncp,dx_ncp));
                auto alpha_tr = zeta*delta/norm_dxncp;
                auto alpha = std::min(alpha_tr,alpha_x_qn);
                if(alpha < Real(1.0)) {
                    X::scal(alpha,dx_ncp);
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
                std::tie(augsys_qn_err,augsys_qn_iter) =
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
                augsys_qn_iter_total+=augsys_qn_iter;
                augsys_iter_total+=augsys_qn_iter;

                // Find the Newton shift, dx_dnewton = dx_newton-dx_ncp
                X_Vector & dx_dnewton = x0.first;

                // Find the Newton step
                X::copy(dx_ncp,dx_n);
                X::axpy(Real(1.),dx_dnewton,dx_n);

                // Safeguard the Newton step 
                alpha_x_qn=std::min(Real(1.),safeguard(dx_ncp,dx_dnewton,zeta));

                // If the dx_n is smaller than zeta delta and we don't have to
                // safeguard, then return it as the quasi-normal step
                Real norm_dxn = sqrt(X::innr(dx_n,dx_n));
                if(norm_dxn <= zeta*delta && alpha_x_qn >= Real(1.))return;

                // Otherwise, compute the dogleg step.  In order to accomplish
                // this, we need to find theta so that
                //
                // || dx_ncp + theta dx_dnewton || = zeta*delta
                //
                // and then set
                //
                // dx_n = dx_ncp + min(theta,alpha_x_qn) dx_dnewton.
                Real aa = X::innr(dx_dnewton,dx_dnewton);
                Real bb = Real(2.) * X::innr(dx_dnewton,dx_ncp);
                Real cc = norm_dxncp*norm_dxncp - zeta*zeta*delta*delta;
                auto roots = quad_equation(aa,bb,cc);
                Real theta = roots[0] > roots[1] ? roots[0] : roots[1];
                alpha = std::min(theta,alpha_x_qn);
                X::copy(dx_ncp,dx_n);
                X::axpy(alpha,dx_dnewton,dx_n);
            }
            
            // Sets the tolerances for projecting 
            //
            // grad f(x) + g'(x)*y + H dx_n
            //
            // into the null space of g'(x).
            struct NullspaceProjForGradLagPlusHdxnManipulator
                : public GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public:
                // Disallow constructors
                NO_DEFAULT_COPY_ASSIGNMENT(
                    NullspaceProjForGradLagPlusHdxnManipulator)

                // Grab the states and fns on construction
                NullspaceProjForGradLagPlusHdxnManipulator(
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}

                // Evalulate the manipulator
                void eval(
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_pg = state.xi_pg;
                    Real const & delta = state.delta;
                    Real & augsys_pg_err_target = state.augsys_pg_err_target;

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

                    // Save this desired error
                    augsys_pg_err_target=eps;
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
                Real & augsys_pg_err = state.augsys_pg_err;
                Natural & augsys_pg_iter = state.augsys_pg_iter;
                Natural & augsys_pg_iter_total = state.augsys_pg_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;

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
                std::tie(augsys_pg_err,augsys_pg_iter) =
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
                augsys_pg_iter_total+=augsys_pg_iter;
                augsys_iter_total+=augsys_pg_iter;

                // Copy out the solution
                X::copy(x0.first,W_gradpHdxn);
            }
            
            // Sets the tolerances for the nullspace projector that projects
            // the current direction in the projected Krylov method. 
            struct NullspaceProjForKrylovMethodManipulator
                : public GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public:
                // Disallow constructors
                NO_DEFAULT_COPY_ASSIGNMENT(
                    NullspaceProjForKrylovMethodManipulator)

                // Grab the states and fns on construction
                explicit NullspaceProjForKrylovMethodManipulator (
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}

                // Evalulate the manipulator
                void eval(
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_proj = state.xi_proj;
                    Real& augsys_proj_err_target=state.augsys_proj_err_target;

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

                    // Save this desired error
                    augsys_proj_err_target=eps;
                }
            };
            
            // Nullspace projector that projects the current direction in the
            // projected Krylov method. 
            struct NullspaceProjForKrylovMethod: public Operator <Real,XX,XX> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public:
                NullspaceProjForKrylovMethod(
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}
               
                // Project dx_t into the nullspace of g'(x)
                void eval(
                    X_Vector const & dx_t_uncorrected,
                    X_Vector & result
                ) const{
                    // Create some shortcuts
                    X_Vector const & x=state.x;
                    Y_Vector const & y=state.y;
                    Natural const & augsys_iter_max=state.augsys_iter_max;
                    Natural const & augsys_rst_freq=state.augsys_rst_freq;
                    Real & augsys_proj_err = state.augsys_proj_err;
                    Natural & augsys_proj_iter = state.augsys_proj_iter;
                    Natural & augsys_proj_iter_total = 
                        state.augsys_proj_iter_total;
                    Natural & augsys_iter_total = state.augsys_iter_total;

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
                    auto iter = Natural(0);
                    std::tie(augsys_proj_err,iter) =
                        Optizelle::gmres <Real,XXxYY> (
                            AugmentedSystem(state,fns,x),
                            b0,
                            Real(1.), // Overwritten by the manipulator
                            augsys_iter_max,
                            augsys_rst_freq,
                            PAugSys_l,
                            PAugSys_r,
                            NullspaceProjForKrylovMethodManipulator(state,fns),
                            x0 
                        );
                    augsys_proj_iter+=iter;
                    augsys_proj_iter_total+=iter;
                    augsys_iter_total+=iter;

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
                ScalarValuedFunction <Real,XX> const & f=*(fns.f);
                ScalarValuedFunctionModifications <Real,XX> const & f_mod
                    = *(fns.f_mod);
                auto const & safeguard = *(fns.safeguard); 
                X_Vector const & x=state.x;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & W_gradpHdxn=state.W_gradpHdxn;
                Real const & delta = state.delta;
                Real const & eps_krylov=state.eps_krylov;
                Natural const & krylov_iter_max=state.krylov_iter_max;
                Natural const & krylov_orthog_max=state.krylov_orthog_max;
                auto const & failed_safeguard_max = state.failed_safeguard_max;
                X_Vector & dx_t_uncorrected=state.dx_t_uncorrected;
                X_Vector & dx_tcp_uncorrected=state.dx_tcp_uncorrected;
                Real & krylov_rel_err=state.krylov_rel_err;
                Natural & krylov_iter=state.krylov_iter;
                Natural & krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                Natural & augsys_proj_iter = state.augsys_proj_iter;
                auto & failed_safeguard = state.failed_safeguard;
                auto & failed_safeguard_total = state.failed_safeguard_total;
                auto & alpha_x = state.alpha_x;
                
                // Setup the Hessian operator and allocate memory for the
                // Cauchy point.
                typename Unconstrained <Real,XX>::Algorithms::HessianOperator
                    H(fns,x);

                // Find the quantity - W (g + H dxn).  We use this as the
                // RHS in the linear system solve.
                X_Vector minus_W_gradpHdxn(X::init(x));
                    X::copy(W_gradpHdxn,minus_W_gradpHdxn);
                    X::scal(Real(-1.),minus_W_gradpHdxn);

                // Keep track of the residual errors
                Real residual_err0(std::numeric_limits <Real>::quiet_NaN());
                Real residual_err(std::numeric_limits <Real>::quiet_NaN());

                // Create the simplified safeguard function
                auto simplified_safeguard =
                    SafeguardSimplified <Real,XX>(std::bind(
                        safeguard,
                        std::placeholders::_1,
                        std::placeholders::_2,
                        Real(1.)));

                // Make sure to zero out our iteration counter for the
                // nullspace projection.  We'll do several iterations of
                // CD and we accumulate this number as we go
                augsys_proj_iter=0;

                // Find the trial step 
                truncated_cg(
                    H,
                    minus_W_gradpHdxn,
                    NullspaceProjForKrylovMethod(state,fns), // Add in PH?
                    eps_krylov,
                    krylov_iter_max,
                    krylov_orthog_max,
                    delta,
                    dx_n,
                    true,
                    failed_safeguard_max,
                    simplified_safeguard,
                    dx_t_uncorrected,
                    dx_tcp_uncorrected,
                    residual_err0,
                    residual_err,
                    krylov_iter,
                    krylov_stop,
                    failed_safeguard,
                    alpha_x);

                // Calculate the Krylov error
                krylov_rel_err = residual_err / residual_err0;
                krylov_iter_total += krylov_iter;

                // Keep track of the number of failed safeguard steps
                failed_safeguard_total+=failed_safeguard;
            }
            
            // Sets the tolerances for the computation of the tangential
            // step.
            struct TangentialStepManipulator
                : public GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public: 
                // Disallow constructors
                NO_DEFAULT_COPY_ASSIGNMENT(TangentialStepManipulator)

                // Grab the states and fns on construction
                TangentialStepManipulator (
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : state(state_), fns(fns_) {}

                // Evalulate the manipulator
                void eval(
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    X_Vector const & dx_n=state.dx_n;
                    Real const & xi_tang = state.xi_tang;
                    Real const & delta = state.delta;
                    Real& augsys_tang_err_target=state.augsys_tang_err_target;

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

                    // Save this desired error
                    augsys_tang_err_target=eps;
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
                Real & augsys_tang_err = state.augsys_tang_err;
                Natural & augsys_tang_iter = state.augsys_tang_iter;
                Natural & augsys_tang_iter_total = state.augsys_tang_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;

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
                std::tie(augsys_tang_err,augsys_tang_iter) =
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
                augsys_tang_iter_total+=augsys_tang_iter;
                augsys_iter_total+=augsys_tang_iter;

                // Copy out the tangential step
                X::copy(x0.first,dx_t);
            }
            
            // Sets the tolerances for the computation of the equality 
            // multiplier.
            struct EqualityMultiplierStepManipulator
                : public GMRESManipulator <Real,XXxYY> {
            private:
                typename State::t & state;
                typename Functions::t const & fns;
            public:
                // Disallow constructors
                NO_DEFAULT_COPY_ASSIGNMENT(EqualityMultiplierStepManipulator)

                // Grab the states and fns on construction
                EqualityMultiplierStepManipulator(
                    typename State::t & state_,
                    typename Functions::t const & fns_
                ) : GMRESManipulator <Real,XXxYY>(), state(state_), fns(fns_) {}

                // Evalulate the manipulator
                void eval(
                    Natural const & iter,
                    typename XXxYY <Real>::Vector const & xx,
                    typename XXxYY <Real>::Vector const & bb,
                    Real & eps
                ) const {
                    // Create some shortcuts
                    Real const & xi_lmh = state.xi_lmh;
                    Real const & xi_lmg = state.xi_lmg;
                    Real & augsys_lmh_err_target= state.augsys_lmh_err_target;
                
                    // Find the norm of the gradient of the Lagrangian.
                    // Sometimes, this is -grad L(x+dx,y).  Sometimes, this
                    // is -grad L(x,y).  In both cases, we just look at the
                    // first element of the RHS.
                    Real norm_grad = sqrt(X::innr(bb.first,bb.first));

                    // The bound is
                    // min( xi_lmg, xi_lmh || grad f(x) + g'(x)*y ||)
                    eps = xi_lmg < norm_grad*xi_lmh ? xi_lmg : norm_grad*xi_lmh;

                    // Save this desired error
                    augsys_lmh_err_target=eps;
                }
            };

            // Finds the equality multiplier at the current iterate 
            static void findEqualityMultiplier(
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
                Real & augsys_lmh_err = state.augsys_lmh_err;
                Natural & augsys_lmh_iter = state.augsys_lmh_iter;
                Natural & augsys_lmh_iter_total = state.augsys_lmh_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;

                // Find the gradient modifications for the equality multiplier
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

                // Solve the augmented system for the initial equality 
                // multiplier 
                std::tie(augsys_lmh_err,augsys_lmh_iter) =
                    Optizelle::gmres <Real,XXxYY> (
                        AugmentedSystem(state,fns,x),
                        b0,
                        Real(1.), // This will be overwritten by the manipulator
                        augsys_iter_max,
                        augsys_rst_freq,
                        PAugSys_l,
                        PAugSys_r,
                        EqualityMultiplierStepManipulator(state,fns),
                        x0 
                    );
                augsys_lmh_iter_total+=augsys_lmh_iter;
                augsys_iter_total+=augsys_lmh_iter;

                // Find the equality multiplier based on this step
                Y::axpy(Real(1.),x0.second,y);
            }
            
            // Finds the equality multiplier step 
            static void findEqualityMultiplierStep(
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
                Real & augsys_lmh_err = state.augsys_lmh_err;
                Natural & augsys_lmh_iter = state.augsys_lmh_iter;
                Natural & augsys_lmh_iter_total = state.augsys_lmh_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;

                // x_p_dx <- x + dx
                X_Vector x_p_dx(X::init(x));
                    X::copy(x,x_p_dx);
                    X::axpy(Real(1.),dx,x_p_dx);

                // grad_xpdx <- L(x+dx,y) = grad f(x+dx) + g'(x+dx)*y
                X_Vector grad_xpdx(X::init(x));
                    f.grad(x_p_dx,grad_xpdx);

                // Find the gradient modifications for the equality multiplier
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

                // Solve the augmented system for the equality multiplier step 
                std::tie(augsys_lmh_err,augsys_lmh_iter) =
                    Optizelle::gmres <Real,XXxYY> (
                        AugmentedSystem(state,fns,x),
                        b0,
                        Real(1.), // This will be overwritten by the manipulator
                        augsys_iter_max,
                        augsys_rst_freq,
                        PAugSys_l,
                        PAugSys_r,
                        EqualityMultiplierStepManipulator(state,fns),
                        x0 
                    );
                augsys_lmh_iter_total+=augsys_lmh_iter;
                augsys_iter_total+=augsys_lmh_iter;

                // Restore our current iterate
                X::copy(x_save,x);

                // Copy out the equality multiplier step
                Y::copy(x0.second,dy);
            }
            
            // Does a check on how far off the equality multiplier is 
            static Real equalityMultiplierCheck(
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
                Real & augsys_lmh_err = state.augsys_lmh_err;
                Natural & augsys_lmh_iter = state.augsys_lmh_iter;
                Natural & augsys_lmh_iter_total = state.augsys_lmh_iter_total;
                Natural & augsys_iter_total = state.augsys_iter_total;

                // grad_x <- L(x,y) = grad f(x) + g'(x)*y
                X_Vector grad_x(X::init(x));
                    f.grad(x,grad_x);

                // Find the gradient modifications for the equality multiplier
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

                // Solve the augmented system for the equality multiplier step 
                std::tie(augsys_lmh_err,augsys_lmh_iter) =
                    Optizelle::gmres <Real,XXxYY> (
                        AugmentedSystem(state,fns,x),
                        b0,
                        Real(1.), // This will be overwritten by the manipulator
                        augsys_iter_max,
                        augsys_rst_freq,
                        PAugSys_l,
                        PAugSys_r,
                        EqualityMultiplierStepManipulator(state,fns),
                        x0 
                    );
                augsys_lmh_iter_total+=augsys_lmh_iter;
                augsys_iter_total+=augsys_lmh_iter;

                // Copy out the equality multiplier step
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

                // If the predicted reduction is small or negative, update the
                // penalty parameter.  Essentially, we're going to force the
                // predicted reduction to be positive as long as the
                // quasinormal step improved feasibility and the tangential
                // step obtained reduction.  Our choice here guarantees
                //
                // pred >= (rho/2) (|| g(x) || - || g'(x) dx_n + g(x) ||)
                //
                // Further, when we update rho_new > rho_old thanks to rho_bar.
                if( pred < (rho_old/Real(2.))
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

                // Save the old equality multiplier
                Y_Vector y_old(Y::init(y));
                    Y::copy(y,y_old);

                // Determine y + dy
                Y::axpy(Real(1.),dy,y);

                // Determine the merit function at x and x+dx
                Real merit_x = f_mod.merit(x,f_x);
                f_xpdx = f.eval(x_p_dx);
                Real merit_xpdx = f_mod.merit(x_p_dx,f_xpdx);

                // Restore the old equality multiplier
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
                auto const & safeguard = *(fns.safeguard);
                auto const & gradmod=*(fns.gradmod);
                auto const & absrel = *(fns.absrel);
                X_Vector const & x=state.x;
                X_Vector const & dx_n=state.dx_n;
                X_Vector const & dx_tcp_uncorrected
                    =state.dx_tcp_uncorrected;
                Real const & xi_4=state.xi_4;
                Real const & eta0=state.eta0;
                Real const & eps_dx=state.eps_dx;
                Real const & norm_dxtyp=state.norm_dxtyp;
                Real const & rho_old=state.rho_old;
                Natural const & history_reset=state.history_reset;
                auto const & eps_constr = state.eps_constr;
                auto const & norm_gxtyp = state.norm_gxtyp;
                auto const & W_gradpHdxn = state.W_gradpHdxn;
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
                auto & alpha_x = state.alpha_x;
                auto & dx_t=state.dx_t;

                // Create a single temporary vector
                X_Vector x_tmp1(X::init(x));

                // Continue to look for a step until our actual vs. predicted
                // reduction is good.
                rejected_trustregion=0;
                while(true) {
                    // Manipulate the state if required
                    smanip.eval(fns,state,OptimizationLocation::BeforeGetStep);

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
                            
                            // See if we need to modify the gradient.  If we do
                            // recompute.
                            auto norm_gx=sqrt(Y::innr(g_x,g_x)); 
                            auto gx_reduction =
                                log10(absrel(norm_gxtyp))-log10(norm_gx);
                            auto gx_converged =
                                norm_gx < eps_constr * absrel(norm_gxtyp);
                            if(gradmod(W_gradpHdxn,gx_reduction,gx_converged))
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

                            // Safeguard the tangential step.  In theory, this
                            // step should already be feasible, but with
                            // inexactness, dx_t and dx_t_uncorrected may be
                            // moderately different, so we need to safeguard
                            // again.
                            auto alpha_safeguard = std::min(Real(1.),
                                safeguard(dx_n,dx_t,Real(1.)));
                            if(alpha_safeguard < Real(1.)) {
                                alpha_x *= alpha_safeguard;
                                X::scal(alpha_safeguard,dx_t);
                            }

                            // Find the primal step
                            X::copy(dx_n,dx);
                            X::axpy(Real(1.),dx_t,dx);

                            // Find g'(x)dx_t
                            g.p(x,dx_t,gpxdxt);

                            // Find the equality multiplier step
                            findEqualityMultiplierStep(fns,state);

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
                                break;
                            }
                        }
                    }

                    // Manipulate the state if required
                    smanip.eval(fns,state,
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
                    smanip.eval(fns,state,
                        OptimizationLocation::AfterRejectedTrustRegion);

                    // Alternatively, check if the step becomes so small
                    // that we're not making progress or if we have a step
                    // with NaNs in it.  In this case, break and allow the
                    // stopping conditions to terminate optimization.  We use a
                    // zero length step so that we do not modify the current
                    // iterate.
                    Real norm_dx = sqrt(X::innr(dx,dx));
                    if( norm_dx < eps_dx * absrel(norm_dxtyp) ||
                        norm_dx != norm_dx
                    ) {
                        X::zero(dx);
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
                auto const & absrel = *(fns.absrel);
                Y_Vector const & g_x = state.g_x;
                Real const & eps_constr=state.eps_constr;
                Real const & norm_gxtyp=state.norm_gxtyp; 
                StoppingCondition::t & opt_stop=state.opt_stop;
                
                // Prevent convergence unless the infeasibility is small. 
                Real norm_gx=sqrt(Y::innr(g_x,g_x));
                if( opt_stop==StoppingCondition::GradientSmall &&
                    !(norm_gx < eps_constr * absrel(norm_gxtyp)) 
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
                void eval(
                    typename ProblemClass::Functions::t const & fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t const & loc
                ) const {
                    // Call the user define manipulator
                    smanip.eval(fns_,state_,loc);

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
                    auto const & multsolve=*(fns.multsolve);
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
                        g.eval(x,g_x);
                        norm_gxtyp = sqrt(Y::innr(g_x,g_x));

                        // Find the initial equality multiplier and then update
                        // the gradient and merit function.
                        findEqualityMultiplier(fns,state);

                        // In addition, update the norm of gradient and
                        // typical gradient since we've modified the equality 
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

                    case OptimizationLocation::AfterStepBeforeGradient:
                        // Make sure we update our cached value of g(x) 
                        g.eval(x,g_x);
                        break;

                    case OptimizationLocation::AfterGradient:
                        // If we're also computing with a primal-dual interior
                        // point method, we likely would have modified our
                        // inequality multiplier prior to the gradient, which
                        // necessitates a new equality multiplier computation.
                        if(multsolve())
                            findEqualityMultiplier(fns,state);
                        break;

                    case OptimizationLocation::AfterCheckStop:
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
                EmptyManipulator <EqualityConstrained <Real,XX,YY> > smanip;

                // Minimize the problem
                getMin(msg,fns,state,smanip);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state,
                StateManipulator <EqualityConstrained <Real,XX,YY> > const &
                    smanip
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
        NO_CONSTRUCTORS(InequalityConstrained)

        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
            // Disallow constructors
            NO_CONSTRUCTORS(State)

            // Internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                NO_DEFAULT_COPY_ASSIGNMENT(t)

                // Inequality multiplier (dual variable or Lagrange multiplier)
                Z_Vector z;
                
                // Step in the inequality multiplier 
                Z_Vector dz;

                // The inequality constraint evaluated at x.  In theory,
                // we can always just evaluate this when we need it.  However,
                // we require its computation both in the gradient as well as
                // Hessian calculations.  More specifically, when computing
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

                // Rate that we decrease the interior point parameter
                Real sigma;

                // How close we move to the boundary during a single step
                Real gamma;

                // Amount we truncate dz in order to maintain feasibility
                // of the inequality multiplier
                Real alpha_z;

                // Function diagnostics on h
                FunctionDiagnostics::t h_diag;

                // Vector space diagnostics on Z
                VectorSpaceDiagnostics::t z_diag;

                // Initialization constructors
                t(X_Vector const & x_user,Z_Vector const & z_user) :
                    Unconstrained <Real,XX>::State::t(x_user),
                        //---z0---
                        // mu inv(L(h(x))) e
                        //---z1---
                    z(Z::init(z_user)),
                    dz(
                        //---dz0---
                        Z::init(z_user)
                        //---dz1---
                    ),
                    h_x(
                        //---h_x0---
                        Z::init(z_user)
                        //---h_x1---
                    ),
                    mu(
                        //---mu0---
                        1.0
                        //---mu1---
                    ),
                    mu_est(
                        //---mu_est0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---mu_est1---
                    ),
                    mu_typ(
                        //---mu_typ0---
                        std::numeric_limits<Real>::quiet_NaN()
                        //---mu_typ1---
                    ),
                    eps_mu(
                        //---eps_mu0---
                        1e-8
                        //---eps_mu1---
                    ),
                    sigma(
                        //---sigma0---
                        0.1
                        //---sigma1---
                    ),
                    gamma(
                        //---gamma0---
                        0.99
                        //---gamma1---
                    ),
                    alpha_z(
                        //---alpha_z0---
                        std::numeric_limits <Real>::quiet_NaN() 
                        //---alpha_z1---
                    ),
                    h_diag(
                        //---h_diag0---
                        FunctionDiagnostics::NoDiagnostics
                        //---h_diag1---
                    ),
                    z_diag(
                        //---z_diag0---
                        VectorSpaceDiagnostics::NoDiagnostics
                        //---z_diag1---
                    )
                {
                        Z::copy(z_user,z);
                }
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };
                
            // Check that we have a valid set of parameters.  
            static void check_(Messaging const & msg,t const & state) {
                // Use this to build an error message
                std::stringstream ss;
                    
                    //---z_valid0---
                    // Any
                    //---z_valid1---
                    
                    //---dz_valid0---
                    // Any
                    //---dz_valid1---
                    
                    //---h_x_valid0---
                    // Any
                    //---h_x_valid1---
                
                // Check that the interior point parameter is positive
                if(!(
                    //---mu_valid0---
                    state.mu > Real(0.)
                    //---mu_valid1---
                )) 
                    ss << "The interior point parameter must be positive: " 
                        "mu = " << state.mu;
                
                // Check that the estimated interior point parameter is a 
                // non-nan past iteration 1.
                else if(!(
                    //---mu_est_valid0---
                    state.mu_est == state.mu_est || state.iter == 1
                    //---mu_est_valid1---
                )) 
                    ss << "The estimated interior point parameter must be "
                        "number: mu_est = " << state.mu_est;

                // Check that the typical interior point parameter is positive 
                // or if we're on the first iteration, we allow NaN.
                else if(!(
                    //---mu_typ_valid0---
                    state.mu_typ > Real(0.) || state.iter==1 
                    //---mu_typ_valid1---
                ))
                    ss << "The typical interior point parameter must be "
                        "positive:  mu_typ = " << state.mu_typ;

                // Check that the interior point stopping tolerance is positive 
                else if(!(
                    //---eps_mu_valid0---
                    state.eps_mu > Real(0.)
                    //---eps_mu_valid1---
                ))
                    ss << "The interior point stopping tolerance must be "
                        "positive: eps_mu = " << state.eps_mu;

                // Check that the reduction in the interior point parameter
                // is between 0 and 1.
                else if(!(
                    //---sigma_valid0---
                    state.sigma > Real(0.) && state.sigma < Real(1.)
                    //---sigma_valid1---
                ))
                    ss << "The reduction in the interior point parameter "
                        "must be between 0 and 1: sigma = " << state.sigma;

                // Check that the fraction to the boundary is between 0 and 1. 
                else if(!(
                    //---gamma_valid0---
                    state.gamma > Real(0.) && state.gamma < Real(1.)
                    //---gamma_valid1---
                ))
                    ss << "The fraction to the boundary must be between " 
                        "0 and 1: gamma= " << state.gamma;
                    
                    //---alpha_z_valid0---
                    // Any
                    //---alpha_z_valid1---
                    
                    //---h_diag_valid0---
                    // Any
                    //---h_diag_valid1---
                    
                    //---z_diag_valid0---
                    // Any
                    //---z_diag_valid1---

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
            NO_CONSTRUCTORS(Restart)
       
            // Create some type shortcuts 
            typedef typename RestartPackage <Real>::t Reals;
            typedef typename RestartPackage <Natural>::t Naturals;
            typedef typename RestartPackage <std::string>::t Params;
            typedef typename RestartPackage <X_Vector>::t X_Vectors;
            typedef typename RestartPackage <Z_Vector>::t Z_Vectors;
            
            // Checks whether we have a valid real 
            static bool is_real(
                typename RestartPackage <Real>::tuple const & item
            ) {
                if( Unconstrained <Real,XX>::Restart::is_real(item) ||
                    item.first == "mu" ||
                    item.first == "mu_est" ||
                    item.first == "mu_typ" ||
                    item.first == "eps_mu" ||
                    item.first == "sigma" ||
                    item.first == "gamma" ||
                    item.first == "alpha_z"
                )
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid natural number
            static bool is_nat(
                typename RestartPackage <Natural>::tuple const & item
            ) {
                if( Unconstrained <Real,XX>::Restart::is_nat(item))
                    return true;
                else
                    return false;
            }
           
            // Checks whether we have a valid parameter 
            static bool is_param(
                typename RestartPackage <std::string>::tuple const & item
            ){
                if( Unconstrained <Real,XX>::Restart::is_param(item) ||
                    (item.first=="h_diag" &&
                        FunctionDiagnostics::is_valid(item.second)) ||
                    (item.first=="z_diag" &&
                        VectorSpaceDiagnostics::is_valid(item.second)) 
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid variable
            static bool is_x(
                typename RestartPackage <X_Vector>::tuple const & item
            ) {
                if( Unconstrained <Real,XX>::Restart::is_x(item))
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid inequality multiplier 
            static bool is_z(
                typename RestartPackage <Z_Vector>::tuple const & item
            ) {
                if( item.first == "z" ||
                    item.first == "dz" ||
                    item.first == "h_x"
                )
                    return true;
                else
                    return false;
            }

            // Checks whether we have valid labels
            static void checkItems(
                Messaging const & msg,
                Reals const & reals,
                Naturals const & nats,
                Params const & params,
                X_Vectors const & xs,
                Z_Vectors const & zs
            ) {
                Utility::checkItems <Real> (
                    msg,is_real,reals," real name: ");
                Utility::checkItems <Natural> (
                    msg,is_nat,nats," natural name: ");
                Utility::checkItems <std::string> (
                    msg,is_param,params," paramater: ");
                Utility::checkItems <X_Vector> (
                    msg,is_x,xs," variable name: ");
                Utility::checkItems <Z_Vector> (
                    msg,is_z,zs," inequality multiplier name: ");
            }
            
            // Copy out the inequality multipliers 
            static void stateToVectors(
                typename State::t & state, 
                X_Vectors & xs,
                Z_Vectors & zs
            ) {
                zs.emplace_back("z",std::move(state.z));
                zs.emplace_back("dz",std::move(state.dz));
                zs.emplace_back("h_x",std::move(state.h_x));
            }
            
            // Copy out the scalar information
            static void stateToScalars(
                typename State::t & state,
                Reals & reals,
                Naturals & nats,
                Params & params
            ) {
                // Copy in all the real numbers
                reals.emplace_back("mu",std::move(state.mu));
                reals.emplace_back("mu_est",std::move(state.mu_est));
                reals.emplace_back("mu_typ",std::move(state.mu_typ));
                reals.emplace_back("eps_mu",std::move(state.eps_mu));
                reals.emplace_back("sigma",std::move(state.sigma));
                reals.emplace_back("gamma",std::move(state.gamma));
                reals.emplace_back("alpha_z",std::move(state.alpha_z));

                // Copy in all of the parameters
                params.emplace_back("h_diag",
                    FunctionDiagnostics::to_string(state.h_diag));
                params.emplace_back("z_diag",
                    VectorSpaceDiagnostics::to_string(state.z_diag));
            }
            
            // Copy in inequality multipliers 
            static void vectorsToState(
                typename State::t & state,
                X_Vectors & xs,
                Z_Vectors & zs
            ) {
                for(typename Z_Vectors::iterator item = zs.begin();
                    item!=zs.end();
                    item++
                ){
                    if(item->first=="z")
                        state.z = std::move(item->second);
                    else if(item->first=="dz")
                        state.dz = std::move(item->second);
                    else if(item->first=="h_x")
                        state.h_x = std::move(item->second);
                }
            }
            
            // Copy in the scalar information
            static void scalarsToState(
                typename State::t & state,
                Reals & reals,
                Naturals & nats,
                Params & params
            ) { 
                // Copy in any reals 
                for(typename Reals::iterator item = reals.begin();
                    item!=reals.end();
                    item++
                ){
                    if(item->first=="mu")
                        state.mu=std::move(item->second);
                    else if(item->first=="mu_est")
                        state.mu_est=std::move(item->second);
                    else if(item->first=="mu_typ")
                        state.mu_typ=std::move(item->second);
                    else if(item->first=="eps_mu")
                        state.eps_mu=std::move(item->second);
                    else if(item->first=="sigma")
                        state.sigma=std::move(item->second);
                    else if(item->first=="gamma")
                        state.gamma=std::move(item->second);
                    else if(item->first=="alpha_z")
                        state.alpha_z=std::move(item->second);
                } 
                    
                // Next, copy in any parameters 
                for(typename Params::iterator item = params.begin();
                    item!=params.end();
                    item++
                ){
                    if(item->first=="h_diag")
                        state.h_diag=FunctionDiagnostics::from_string(
                            item->second);
                    else if(item->first=="z_diag")
                        state.z_diag=VectorSpaceDiagnostics::from_string(
                            item->second);
                }
            }

            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Z_Vectors & zs,
                Reals & reals,
                Naturals & nats,
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
                Naturals & nats,
                Params & params
            ) {
                // Check the user input 
                checkItems(msg,reals,nats,params,xs,zs);

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
            NO_CONSTRUCTORS(Functions)

            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t)

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
                NO_COPY_ASSIGNMENT(InequalityModifications)

            private:
                // Underlying modification.  This takes control of the memory
                std::unique_ptr <
                    Optizelle::ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Inequality constraint.
                Optizelle::VectorValuedFunction <Real,XX,ZZ> const & h;
                
                // Inequality multiplier
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
                    typename Functions::t const & fns,
                    typename State::t const & state,
                    std::unique_ptr <
                        ScalarValuedFunctionModifications<Real,XX>
                    > && f_mod_
                ) : f_mod(std::move(f_mod_)),
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
                        h.eval(x,hx_merit);
                        
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

            // Safeguard search for inequality constraints
            //
            // max { h'(x) alpha dx_dir +
            //              zeta gamma h(x) + h'(x) dx_base >=0 }
            //
            //    or srch(h'(x) dx_dir, zeta gamma h(x) + h'(x) dx_base)
            //
            // Basically, the first safeguard is the fraction to the boundary
            // rule that insures
            //
            // h(x + dx_base + alpha dx_dir)>= (1-gamma zeta) h(x).
            //
            // When we're not doing a composite-step method, we just set zeta
            // to 1.  Really, we need it to preserve some slack between the
            // quasinormal and tangential steps.  For purely inequality
            // constrained problems, we push as far as possible each time.
            static Real inequalitySafeguard(
                typename Functions::t const & fns,
                typename State::t const & state,
                X_Vector const & dx_base,
                X_Vector const & dx_dir,
                Real const & zeta
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                auto const & x = state.x;
                auto const & z = state.z;
                auto const & h_x = state.h_x;
                auto const & gamma = state.gamma;

                // hp_dxdir <- h'(x)dx_dir
                auto hp_dxdir = Z::init(z);
                h.p(x,dx_dir,hp_dxdir);

                // hp_dxbase <- h'(x)dx_base
                auto hp_dxbase = Z::init(z);
                h.p(x,dx_base,hp_dxbase);

                // base <- gamma zeta h(x) + h'(x)dx_base
                auto base = Z::init(z);
                Z::copy(hp_dxbase,base); 
                Z::axpy(gamma*zeta,h_x,base);

                // alpha_safeguard =
                // max{ h'(x) alpha dx_dir + gamma zeta h(x) + h'(x)dx_base >=0}
                auto alpha_safeguard = Z::srch(hp_dxdir,base); 

                // Return the amount we need to safeguard
                return alpha_safeguard;
            }
           
            // Gradient step modification for the inequality constraints.
            // Basically, we can run into trouble when the gradient passed into
            // the step calculation becomes really small prior to convergence.
            // Essentially, the algorithm thinks that it converged, produces a
            // small step, and sometimes this step is small enough that it
            // exceeds our floating point error and we get bad things like
            // zero, or slightly negative, predicted or actual reduction.  To
            // combat this, we can monitor how small the norm of grad_step
            // really is.  If it's small, then it's likely time that we modify
            // our interior point parameter.  For the equality constrained
            // problem, this means that we need to recompute our projected
            // gradient, which is why we do this here as opposed to the state
            // manipulator.
            static Real inequalityGradStepModification(
                typename Functions::t const & fns,
                typename State::t & state,
                X_Vector const & grad_step,
                Real const & gx_reduction,
                bool const & gx_converged
            ) {
                // Create some shortcuts
                auto const & f_mod = *(fns.f_mod);
                auto const & absrel = *(fns.absrel);
                auto const & x=state.x;
                auto const & grad=state.grad;
                auto const & norm_gradtyp=state.norm_gradtyp;
                auto const & eps_grad=state.eps_grad;
                auto const & mu_typ=state.mu_typ;
                auto const & eps_mu=state.eps_mu;
                auto const & mu_est=state.mu_est;
                auto const & sigma=state.sigma;
                auto const & iter = state.iter;
                auto & mu=state.mu;
                auto & z=state.z;
                auto & dz=state.dz;

                // If mu satisfies the stopping criteria, stop trying to
                // reduce the interior point parameter
                if(std::fabs(mu-absrel(mu_typ)*eps_mu)<absrel(mu_typ)*eps_mu)
                    return false;

                // Find || grad_step ||
                auto norm_grad_step = std::sqrt(X::innr(grad_step,grad_step));
                       
                // Find || grad_stop ||
                auto grad_stop = X::init(x);
                f_mod.grad_stop(x,grad,grad_stop);
                auto norm_grad_stop = std::sqrt(X::innr(grad_stop,grad_stop));

                // We reduce mu when the following has occured
                //
                // 1. We converge grad either locally or globally
                //    (need only one)
                //
                //    a. || grad_step || compared to || grad_typ || is smaller
                //       than mu_typ compared to mu_est
                //
                //    b. || grad_step || < eps_grad || grad_typ ||
                //       
                //    c. || grad_stop || compared to || grad_typ || is smaller
                //       than mu_typ compared to mu_est
                //
                //    d. || grad_stop || < eps_grad || grad_typ ||
                //
                // 2. We converge g(x) locally or globally (need only one)
                //
                //    a. || g(x) || compared to || g(x)_typ || is smaller
                //       than mu_typ compared to mu_est
                //
                //    b. || g(x) || < eps_constr || g(x)_typ ||
                //
                // 3. We converge mu_est locally
                //
                //    |mu - mu_est | < mu
                //
                // 4. We have *not* converged mu globally
                //
                //    |mu - eps_mu mu_typ| >= eps_mu mu_typ
                //
                // 5. We're not on the first iteration.  We really haven't
                //    done anything yet, so the reductions are basically zero
                //    except that the grad_step_reduction may be slight because
                //    of how we calculate things.  Really, we're on the first
                //    iteration, give it some time with the specified value.
                auto grad_step_reduction =
                    log10(absrel(norm_gradtyp))-log10(norm_grad_step);
                auto grad_step_converged =
                    norm_grad_step < eps_grad*absrel(norm_gradtyp);

                auto grad_stop_reduction =
                    log10(absrel(norm_gradtyp))-log10(norm_grad_stop);
                auto grad_stop_converged =
                    norm_grad_stop < eps_grad*absrel(norm_gradtyp);

                auto mu_reduction = log10(absrel(mu_typ)) - log10(mu_est);
                auto mu_est_converged = std::fabs(mu-mu_est) < mu;
                auto mu_converged =
                    std::fabs(mu-absrel(mu_typ)*eps_mu) < absrel(mu_typ)*eps_mu;

                if( iter>1 &&
                    ((grad_step_reduction >=mu_reduction)||grad_step_converged||
                    (grad_stop_reduction >=mu_reduction)||grad_stop_converged)&&
                    ((gx_reduction >= mu_reduction) || gx_converged) &&
                    mu_est_converged &&
                    !mu_converged
                ) {
                    // Reduce mu;
                    mu *= sigma;
                    
                    // If we're doing a log-barrier method, update our
                    // multiplier
                    if(!Algorithms::usePrimalDual(state))
                        Algorithms::findInequalityMultiplierLogBarrier(state);

                    // Notify that we modified things
                    return true;
                }
                else
                    return false;
            }

            // Notify the equality constrained solve that an additional
            // multiplier solve is required
            static bool inequalityMultiplierSolve() {
                return true;
            }

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
                fns.f_mod.reset(new InequalityModifications(
                    fns,state,std::move(fns.f_mod)));

                // Put in the inequality safeguards and step modification
                fns.safeguard = std::make_unique <Safeguard <Real,XX>>(
                    std::bind(
                        inequalitySafeguard,
                        std::cref(fns),
                        std::cref(state),
                        std::placeholders::_1,
                        std::placeholders::_2,
                        std::placeholders::_3));
                fns.gradmod =
                    std::make_unique <GradStepModification <Real,XX>>(
                        std::bind(
                            inequalityGradStepModification,
                            std::cref(fns),
                            std::ref(state),
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3));
                
                // If we also have equality constraints, we need another
                // multiplier solve 
                fns.multsolve = std::make_unique <MultiplierSolve> (
                    inequalityMultiplierSolve);
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
        struct Diagnostics {
            // Disallow constructors
            NO_CONSTRUCTORS(Diagnostics)

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                // Create some shortcuts
                auto const & msg_level = state.msg_level;
                auto const & algorithm_class = state.algorithm_class;

                // Basic information
                out.emplace_back(Utility::atos("mu_est"));

                // More detailed information
                if(msg_level >= 2) {
                    out.emplace_back(Utility::atos("mu"));
                    out.emplace_back(Utility::atos("alpha_x"));
                    if(Algorithms::usePrimalDual(state))
                        out.emplace_back(Utility::atos("alpha_z"));
                    if(algorithm_class!=AlgorithmClass::LineSearch)
                        out.emplace_back(Utility::atos("failed_safe"));
                }
            }

            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getStateHeader_(
                    state,out);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::getStateHeader_(
                    state,out);
            }

            // Gets the state information for output
            static void getState_(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                auto const & mu=state.mu; 
                auto const & mu_est=state.mu_est; 
                auto const & alpha_x =state.alpha_x;
                auto const & alpha_z =state.alpha_z;
                auto const & failed_safeguard=state.failed_safeguard;
                auto const & msg_level = state.msg_level;
                auto const & algorithm_class = state.algorithm_class;

                // Figure out if we're at the absolute beginning of the
                // optimization.
                bool opt_begin = Utility::is_opt_begin <InequalityConstrained> (
                    state);

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                out.emplace_back(Utility::atos(mu_est));
                
                // More detailed information
                if(msg_level >= 2) {
                    out.emplace_back(Utility::atos(mu));

                    if(!opt_begin) {
                        out.emplace_back(Utility::atos(alpha_x));
                        if(Algorithms::usePrimalDual(state))
                            out.emplace_back(Utility::atos(alpha_z));
                        if(algorithm_class!=AlgorithmClass::LineSearch)
                            out.emplace_back(Utility::atos(failed_safeguard));
                    } else {
                        auto nblank = 1;
                        if(Algorithms::usePrimalDual(state))
                            nblank++;
                        if(algorithm_class!=AlgorithmClass::LineSearch)
                            nblank++;
                        for(Natural i=0;i<nblank;i++)
                            out.emplace_back(Utility::blankSeparator);
                    }
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
                Unconstrained <Real,XX>::Diagnostics
                    ::getState_(fns,state,blank,noiter,out);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
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
                Unconstrained <Real,XX>::Diagnostics::getKrylovHeader_(
                    state,out);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::getKrylovHeader_
                    (state,out);
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
                Unconstrained <Real,XX>::Diagnostics
                    ::getKrylov_(state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::getKrylov_(state,blank,out);
            }
           
            // Runs the specified function diagnostics 
            static void checkFunctions_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                X_Vector const & x=state.x;
                Z_Vector const & z=state.z;
                FunctionDiagnostics::t const & h_diag = state.h_diag;
                
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 
                Z_Vector dz(Z::init(z));
                    Z::rand(dz);

                // Run the diagnostics
                switch(h_diag) {
                    case FunctionDiagnostics::FirstOrder:
                        msg.print("Diagnostics on the function h");
                        Optizelle::Diagnostics::derivativeCheck(
                            msg,h,x,dx,dz,"h");
                        Optizelle::Diagnostics::derivativeAdjointCheck(
                            msg,h,x,dx,dz,"h");
                        msg.print("");
                        break;
                    case FunctionDiagnostics::SecondOrder:
                        msg.print("Diagnostics on the function h");
                        Optizelle::Diagnostics::derivativeCheck(
                            msg,h,x,dx,dz,"h");
                        Optizelle::Diagnostics::derivativeAdjointCheck(
                            msg,h,x,dx,dz,"h");
                        Optizelle::Diagnostics::secondDerivativeCheck(
                            msg,h,x,dx,dz,"h");
                        msg.print("");
                        break;
                    case FunctionDiagnostics::NoDiagnostics:
                        break;
                }
            }
            
            // Runs the specified function diagnostics 
            static void checkFunctions(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkFunctions_(
                    msg,fns,state);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::checkFunctions_(
                    msg,fns,state);
            }
            
            // Sets up the equality Hessian operator 
            struct InequalityHessianOperator : public Operator <Real,XX,XX> {
            private:
                // Store the equality constraints 
                VectorValuedFunction <Real,XX,ZZ> const & h;

                // Current iterate
                X_Vector const & x; 

                // Inequality modifications 
                typename Functions::InequalityModifications h_mod;

            public:
                // Remove some constructors
                NO_COPY_ASSIGNMENT(InequalityHessianOperator);

                // Take in the objective and the base point during construction 
                InequalityHessianOperator(
                    typename Functions::t const & fns,
                    typename State::t const & state
                ) :
                    h(*(fns.h)),
                    x(state.x),
                    h_mod(
                        fns,
                        state,
                        std::unique_ptr <
                            ScalarValuedFunctionModifications <Real,XX>
                        > (new ScalarValuedFunctionModifications <Real,XX> ())
                    )
                {}

                // Basic application
                void eval(X_Vector const & dx,X_Vector & result)
                    const
                {
                    // Grab the zero vector in the X space
                    X_Vector zero(X::init(x));
                    X::zero(zero);

                    // Add the equality constraint's contribution to the
                    // Hessian-vector product
                    h_mod.hessvec_step(x,dx,zero,result);
                }
            };
           
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                X_Vector const & x=state.x;
                FunctionDiagnostics::t const & L_diag=state.L_diag;
                
                // Create some random directions for these tests
                X_Vector dx(X::init(x));
                    X::rand(dx); 
                X_Vector dxx(X::init(x));
                    X::rand(dxx); 

                // Create the inequality Hessian operator
                InequalityHessianOperator L(fns,state);

                // Run the diagnostics
                switch(L_diag) {
                    case FunctionDiagnostics::FirstOrder:
                    case FunctionDiagnostics::SecondOrder:
                        msg.print("Diagnostics on the contribution of h to "
                            "the Lagrangian");
                        Optizelle::Diagnostics::operatorSymmetryCheck <Real,XX>(
                            msg,L,dx,dxx,"h'(x)*(Linv(h(x))(h'(x).z)");
                        msg.print("");
                        break;
                }
            }
            
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkLagrangian_(msg,fns,state);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::checkLagrangian_(msg,fns,state);
            }

            // Runs the specified vector space diagnostics 
            static void checkVectorSpace_(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                // Create some shortcuts
                VectorValuedFunction <Real,XX,ZZ> const & h=*(fns.h);
                VectorSpaceDiagnostics::t const & z_diag=state.z_diag;
                X_Vector const & x=state.x;
                Z_Vector const & z=state.z;
               
                // Create some random directions for these tests
                Z_Vector dz(Z::init(z));
                    Z::rand(dz); 
                Z_Vector dzz(Z::init(z));
                    Z::rand(dzz); 
                Z_Vector dzzz(Z::init(z));
                    Z::rand(dzzz); 
                Z_Vector dzzzz(Z::init(z));
                    Z::rand(dzzzz);

                // Run the diagnostics
                switch(z_diag) {
                    case VectorSpaceDiagnostics::Basic:
                        msg.print("Diagnostics on the vector-space Z");
                        Optizelle::Diagnostics::zero_innr <Real,ZZ> (msg,z,"Z");
                        Optizelle::Diagnostics::copy_axpy_innr <Real,ZZ> (
                            msg,dz,"Z");
                        Optizelle::Diagnostics::copy_scal_innr <Real,ZZ> (
                            msg,dz,"Z");
                        msg.print("");
                        break;
                    case VectorSpaceDiagnostics::EuclideanJordan: {

                        // Evaluate h_x
                        Z_Vector h_x(Z::init(z));
                        h.eval(x,h_x);

                        // Run the diagnostics
                        msg.print("Diagnostics on the vector-space Z");
                        Optizelle::Diagnostics::zero_innr <Real,ZZ> (msg,z,"Z");
                        Optizelle::Diagnostics::copy_axpy_innr <Real,ZZ> (
                            msg,dz,"Z");
                        Optizelle::Diagnostics::copy_scal_innr <Real,ZZ> (
                            msg,dz,"Z");

                        Optizelle::Diagnostics::id_prod <Real,ZZ> (msg,dz,"Z");
                        Optizelle::Diagnostics::prod_linv <Real,ZZ>(
                            msg,dz,dzz,"Z");
                        Optizelle::Diagnostics::id_srch <Real,ZZ> (msg,z,"Z");
                        Optizelle::Diagnostics::linv_id_barr <Real,ZZ> (
                            msg,h_x,dz,"Z");
                        Optizelle::Diagnostics::innr_prod_symm <Real,ZZ> (
                            msg,dz,dzz,dzzz,dzzzz,"Z");
                        msg.print("");
                        break;
                    } 
                }
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms)


            // Finds the new inequality multiplier step
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

            // Find the log-barrier inequality multiplier.  Specifically,
            // we set z to be 
            //
            // mu inv(L(h(x))) e
            //
            // In this way,
            //
            // h(x) o z = mu e
            //
            // mu_est = <h(x),z> / <e,e>
            //        = mu
            static void findInequalityMultiplierLogBarrier(
                typename State::t & state
            ) {
                // Create some shortcuts
                Z_Vector const & h_x=state.h_x;
                Real const & mu=state.mu;
                Z_Vector & z=state.z;

                // z_tmp1 <- e
                auto e = Z::init(z);
                Z::id(e);

                // z <- inv(L(h(x))) e
                Z::linv(h_x,e,z);

                // z <- mu inv(L(h(x))) e
                Z::scal(mu,z);
            }
           
            // Assume that dz has already been calculated.  Truncate the step
            // in order to insure that our fraction to the boundary rule has
            // been satisfied
            //
            // z + alpha_z dz >= (1-gamma) z
            //
            // or that
            //
            // gamma z + alpha_z dz >= 0
            static void adjustInequalityMultiplierStep(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                auto const & z=state.z;
                auto const & gamma=state.gamma;
                auto & dz=state.dz;
                auto & alpha_z=state.alpha_z;

                // Find gamma z
                auto gamma_z = Z::init(z);
                Z::copy(z,gamma_z);
                Z::scal(gamma,gamma_z);

                // See how far we can search in the dz direction
                alpha_z = std::min(Real(1.0),Z::srch(dz,gamma_z));

                // Truncate our step if need be 
                if(alpha_z < Real(1.))
                    Z::scal(alpha_z,dz);
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
           
            // Adjust the stopping conditions unless the criteria below are
            // satisfied.
            static void adjustStoppingConditions(
                typename Functions::t const & fns,
                typename State::t & state
            ) {
                // Create some shortcuts
                auto const & absrel = *(fns.absrel);
                Real const & mu=state.mu;
                Real const & mu_est=state.mu_est;
                Real const & mu_typ=state.mu_typ;
                Real const & eps_mu=state.eps_mu;
                Natural const & iter=state.iter;
                StoppingCondition::t & opt_stop=state.opt_stop;

                // If the estimated interior point paramter is negative, exit
                if(mu_est < Real(0.)) {
                    opt_stop=StoppingCondition::InteriorPointInstability;
                    return;
                }

                // Prevent convergence unless
                //
                // 1.  mu has been reduced to the same order as eps_mu *
                //     mu_normalization
                // 2.  mu_est is on the same order as mu
                auto mu_converged =
                    std::fabs(mu-absrel(mu_typ)*eps_mu) <=absrel(mu_typ)*eps_mu;
                auto mu_est_converged = std::fabs(mu-mu_est) <= mu;
                if( opt_stop==StoppingCondition::GradientSmall &&
                    !(mu_converged && mu_est_converged)
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
                Real & alpha_x=state.alpha_x;
                Real & alpha_z=state.alpha_z;

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
                    h.eval(x_tmp1,z_tmp1);

                // z_tmp1=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                alpha_x=Z::srch(z_tmp1,h_x);

                // Determine how far we can go in the dual variable 

                // Find the largest alpha such that
                // alpha dz + z >=0
                alpha_z=Z::srch(dz,z);

                // Figure out how much to shorten the steps, if at all
                alpha_x = alpha_x*gamma>Real(1.) ? Real(1.) : alpha_x*gamma;
                alpha_z = alpha_z*gamma>Real(1.) ? Real(1.) : alpha_z*gamma;

                // Shorten the inequality multiplier step
                Z::scal(alpha_z,dz);

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(alpha_x,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= alpha_x;
            }

            // Determines if we are currently using a primal-dual method.
            // Basically, whenever we solve some kind of second-order system,
            // like on a Newton method, we do primal dual.  Otherwise, we do
            // log-barrier.
            static bool usePrimalDual(
                typename State::t const & state
            ) {
                // Create some shortcuts 
                auto const & algorithm_class = state.algorithm_class;
                auto const & dir = state.dir;

                // Anytime we have a second-order system, we do primal-dual.
                // The UserDefined case is tricky.  Right now, we use this
                // for the composite-step SQP, so we'll do primal-dual here, but
                // realistically if the user defines a scheme we may want one
                // versus the other.
                if( algorithm_class == AlgorithmClass::TrustRegion ||
                    algorithm_class == AlgorithmClass::UserDefined ||
                    (algorithm_class == AlgorithmClass::LineSearch &&
                    dir == LineSearchDirection::NewtonCG)
                )
                    return true;
                else
                    return false;
            }
            
            // Conduct a line search that preserves positivity of both the
            // primal and dual variables where the primal and dual variables
            // are linked.  Specifically, every iteration, we set
            //
            // dz = -z + inv(L(h(x)))(-h'(x)dx o z + mu e)
            //
            // Therefore, we need
            //
            // z+dz(alpha dx) >= 0
            //
            // Hence, we need
            //
            // inv(L(h(x)))(-h'(x)alpha dx o z + mu e) >= 0
            //
            // Or
            //
            // -alpha h'(x) dx o z + mu e >= 0
            //
            // Simultaneously, we need
            //
            // x + alpha dx >= 0
            //
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
                Real & alpha_x=state.alpha_x;

                // Don't save information into alpha_z since we're ultimately
                // not going to reference it
                Real alpha_z(std::numeric_limits<Real>::quiet_NaN());

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
                h.eval(x_tmp1,z_tmp1);

                // z_tmp2=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                alpha_x=Z::srch(z_tmp1,h_x);

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
                alpha_z=Z::srch(z_tmp2,z_tmp1);

                // Determine the farthest we can go in both variables
                Real alpha0;

                // Only the dual step is restrictive
                if( alpha_x > std::numeric_limits <Real>::max() &&
                    alpha_z <= std::numeric_limits <Real>::max()
                )
                    alpha0 = alpha_z;

                // Only the primal step is restrictive
                else if(alpha_x <= std::numeric_limits <Real>::max() &&
                    alpha_z > std::numeric_limits <Real>::max()
                )
                    alpha0 = alpha_x;

                // Neither step is restrictive
                else if( alpha_x > std::numeric_limits <Real>::max() &&
                    alpha_z > std::numeric_limits <Real>::max()
                )
                    alpha0 = alpha_x;

                // Both steps are restrictive
                else
                    alpha0 = alpha_x < alpha_z ? alpha_x : alpha_z;
                    
                // Next, determine if we need to back off from the
                // boundary or leave the step unchanged.
                alpha0 = alpha0*gamma > Real(1.) ?  Real(1.) : alpha0*gamma;

                // Save the result in alpha_x
                alpha_x = alpha0;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(alpha_x,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= alpha_x;
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
                Real & alpha_x=state.alpha_x;

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
                    h.eval(x_tmp1,z_tmp1);

                // z_tmp1=h(x+dx)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+dx)-h(x)) + h(x) >=0
                alpha_x=Z::srch(z_tmp1,h_x);

                // Figure out how much to shorten the steps, if at all
                alpha_x = alpha_x*gamma>Real(1.) ? Real(1.) : alpha_x*gamma;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(alpha_x,dx);

                // If we're doing a line-search method, make sure
                // we can't line-search past this point
                else
                    state.alpha0 *= alpha_x;
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
                void eval(
                    typename ProblemClass::Functions::t const & fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t const & loc
                ) const {
                    // Call the user define manipulator
                    smanip.eval(fns_,state_,loc);

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
                    Real const & mu = state.mu;
                    Real const & mu_est = state.mu_est;
                    Z_Vector & z=state.z;
                    Z_Vector & h_x=state.h_x;
                    Z_Vector & dz=state.dz;
                    Real & mu_typ = state.mu_typ;

                    switch(loc){
                    case OptimizationLocation::BeforeInitialFuncAndGrad: 

                        // Initialize the value h(x)
                        h.eval(x,h_x);
            
                        // Find the initial inequality multiplier.  Currently,
                        // we choose our initial z to be what it would be
                        // during a log-barrier method. 
                        findInequalityMultiplierLogBarrier(state);

                        // Estimate the interior point parameter
                        estimateInteriorPointParameter(fns,state);

                        // Set the typical value for mu
                        mu_typ=mu_est;
                        break;

                    case OptimizationLocation::BeforeLineSearch:
                    case OptimizationLocation::BeforeActualVersusPredicted:
                        if(usePrimalDual(state)) {
                            // To be sure, I find this a little bit odd.
                            // Basically, the inequality multiplier step dz
                            // depends on dx.  The code appears to work better
                            // when we solve for dz prior to scaling dx from a
                            // line-search when one is required.  For a
                            // trust-region method, computing here doesn't
                            // matter.

                            // Find the inequality multiplier step 
                            findInequalityMultiplierStep(fns,state);
                            
                            // Truncate the inequality multiplier step if need
                            // be 
                            adjustInequalityMultiplierStep(fns,state);
                        }
                        break;

                    case OptimizationLocation::AfterRejectedTrustRegion:
                    case OptimizationLocation::AfterRejectedLineSearch:
                        // After we reject a step, make sure that we take a
                        // zero step in the inequality multiplier.  This is
                        // important in case we exit early due to small steps.
                        Z::zero(dz);
                        break;

                    case OptimizationLocation::BeforeStep:
                        // Take our inequality multiplier step
                        if(usePrimalDual(state))
                            Z::axpy(Real(1.),dz,z);

                        // Find the log-barrier multiplier
                        else
                            findInequalityMultiplierLogBarrier(state);
                        break;

                    case OptimizationLocation::AfterStepBeforeGradient:
                        // Updated our cached copy of h(x)
                        h.eval(x,h_x);

                        // Update the interior point estimate
                        estimateInteriorPointParameter(fns,state);
                        break;

                    // Adjust the interior point parameter and insure that
                    // we do not converge unless the interior point parameter
                    // is small.
                    case OptimizationLocation::AfterCheckStop:
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
                EmptyManipulator <InequalityConstrained <Real,XX,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,fns,state,smanip);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state,
                StateManipulator <InequalityConstrained <Real,XX,ZZ> > const &
                    smanip
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
        NO_CONSTRUCTORS(Constrained)

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
            NO_CONSTRUCTORS(State)
            
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
                NO_DEFAULT_COPY_ASSIGNMENT(t)

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
            NO_CONSTRUCTORS(Restart)
       
            // Create some type shortcuts 
            typedef typename RestartPackage <Real>::t Reals;
            typedef typename RestartPackage <Natural>::t Naturals;
            typedef typename RestartPackage <std::string>::t Params;
            typedef typename RestartPackage <X_Vector>::t X_Vectors;
            typedef typename RestartPackage <Y_Vector>::t Y_Vectors;
            typedef typename RestartPackage <Z_Vector>::t Z_Vectors;
            
            // Checks whether we have a valid real 
            static bool is_real(
                typename RestartPackage <Real>::tuple const & item
            ) {
                if( EqualityConstrained <Real,XX,YY>::Restart::is_real(item) ||
                    InequalityConstrained <Real,XX,ZZ>::Restart::is_real(item)
                )
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid natural number
            static bool is_nat(
                typename RestartPackage <Natural>::tuple const & item
            ) {
                if( EqualityConstrained <Real,XX,YY>::Restart::is_nat(item) ||
                    InequalityConstrained <Real,XX,ZZ>::Restart::is_nat(item)
                )
                    return true;
                else
                    return false;
            }
           
            // Checks whether we have a valid parameter 
            static bool is_param(
                typename RestartPackage <std::string>::tuple const & item
            ){
                if( EqualityConstrained <Real,XX,YY>::Restart::is_param(item) ||
                    InequalityConstrained <Real,XX,ZZ>::Restart::is_param(item)
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid variable
            static bool is_x(
                typename RestartPackage <X_Vector>::tuple const & item
            ) {
                if( EqualityConstrained <Real,XX,YY>::Restart::is_x(item) ||
                    InequalityConstrained <Real,XX,ZZ>::Restart::is_x(item)
                ) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid equality multiplier label
            static bool is_y(
                typename RestartPackage <Y_Vector>::tuple const & item
            ) {
                if( EqualityConstrained <Real,XX,YY>::Restart::is_y(item)) 
                    return true;
                else
                    return false;
            }
            
            // Checks whether we have a valid inequality multiplier 
            static bool is_z(
                typename RestartPackage <Z_Vector>::tuple const & item
            ) {
                if( InequalityConstrained <Real,XX,ZZ>::Restart::is_z(item)) 
                    return true;
                else
                    return false;
            }

            // Checks whether we have valid labels
            static void checkItems(
                Messaging const & msg,
                Reals const & reals,
                Naturals const & nats,
                Params const & params,
                X_Vectors const & xs,
                Y_Vectors const & ys,
                Z_Vectors const & zs
            ) {
                Utility::checkItems <Real> (
                    msg,is_real,reals," real name: ");
                Utility::checkItems <Natural> (
                    msg,is_nat,nats," natural name: ");
                Utility::checkItems <std::string> (
                    msg,is_param,params," paramater: ");
                Utility::checkItems <X_Vector> (
                    msg,is_x,xs," variable name: ");
                Utility::checkItems <Y_Vector> (
                    msg,is_y,ys," equality multiplier name: ");
                Utility::checkItems <Z_Vector> (
                    msg,is_z,zs," inequality multiplier name: ");
            }
            
            // Release the data into structures controlled by the user 
            static void release(
                typename State::t & state,
                X_Vectors & xs,
                Y_Vectors & ys,
                Z_Vectors & zs,
                Reals & reals,
                Naturals & nats,
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
                Naturals & nats,
                Params & params
            ) {
                // Check the user input 
                checkItems(msg,reals,nats,params,xs,ys,zs);

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
            NO_CONSTRUCTORS(Functions)

            // Actual storage of the functions required
            struct t: 
                public EqualityConstrained <Real,XX,YY>::Functions::t,
                public InequalityConstrained <Real,XX,ZZ>::Functions::t
            {
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of unique_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                NO_COPY_ASSIGNMENT(t)

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
        struct Diagnostics {
            // Disallow constructors
            NO_CONSTRUCTORS(Diagnostics)

            // Gets the header for the state information
            static void getStateHeader_(
                typename State::t const & state,
                std::list <std::string> & out
            ) { 
                // Create some shortcuts
                auto const & msg_level = state.msg_level;

                // More detailed information
                if(msg_level >= 2)
                    out.emplace_back(Utility::atos("alpha_x_qn"));
            }
            // Combines all of the state headers
            static void getStateHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getStateHeader_(
                    state,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getStateHeader_(
                    state,out);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::getStateHeader_(
                    state,out);
                Constrained<Real,XX,YY,ZZ>::Diagnostics::getStateHeader_(
                    state,out);
            }

            // Gets the state information for output
            static void getState_(
                typename Functions::t const & fns,
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {

                // Create some shortcuts
                auto const & alpha_x_qn =state.alpha_x_qn;
                auto const & msg_level = state.msg_level;

                // Figure out if we're at the absolute beginning of the
                // optimization.
                bool opt_begin = Utility::is_opt_begin <Constrained> (state);

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // More detailed information
                if(msg_level >= 2) {
                    if(!opt_begin)
                        out.emplace_back(Utility::atos(alpha_x_qn));
                    else {
                        for(Natural i=0;i<1;i++)
                            out.emplace_back(Utility::blankSeparator);
                    }
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
                Unconstrained <Real,XX>::Diagnostics
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics
                    ::getState_(fns,state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::getState_(fns,state,blank,out);
                Constrained <Real,XX,YY,ZZ>::Diagnostics
                    ::getState_(fns,state,blank,out);
            }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                typename State::t const & state,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getKrylovHeader_(
                    state,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getKrylovHeader_(
                    state,out);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::getKrylovHeader_
                    (state,out);
            }

            // Combines all of the Krylov information
            static void getKrylov(
                typename State::t const & state,
                bool const & blank,
                std::list <std::string> & out
            ) {
                Unconstrained <Real,XX>::Diagnostics::getKrylov_(
                    state,blank,out);
                EqualityConstrained <Real,XX,YY>::Diagnostics::getKrylov_(
                    state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics::getKrylov_(
                    state,blank,out);
            }
            
            // Runs the specified function diagnostics 
            static void checkFunctions(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics::checkFunctions_(
                    msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics::checkFunctions_(
                    msg,fns,state);
                InequalityConstrained<Real,XX,ZZ>::Diagnostics::checkFunctions_(
                    msg,fns,state);
            }
            
            // Runs the specified Lagrangian diagnostics 
            static void checkLagrangian(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkLagrangian_(msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics
                    ::checkLagrangian_(msg,fns,state);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::checkLagrangian_(msg,fns,state);
            }
            
            // Runs the specified vector space diagnostics 
            static void checkVectorSpace(
                Messaging const & msg,
                typename Functions::t const & fns,
                typename State::t const & state
            ) {
                Unconstrained <Real,XX>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
                EqualityConstrained <Real,XX,YY>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
                InequalityConstrained <Real,XX,ZZ>::Diagnostics
                    ::checkVectorSpace_(msg,fns,state);
            }
        };
        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
            // Disallow constructors
            NO_CONSTRUCTORS(Algorithms)

            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state
            ){
                // Create an empty state manipulator
                EmptyManipulator <Constrained <Real,XX,YY,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,fns,state,smanip);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                Messaging const & msg,
                typename Functions::t & fns,
                typename State::t & state,
                StateManipulator <Constrained <Real,XX,YY,ZZ> > const & smanip
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
//---Optizelle2---
}
//---Optizelle3---
#endif
