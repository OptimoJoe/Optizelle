#ifndef PEOPT_H
#define PEOPT_H

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
#include<peopt/linalg.h>

namespace peopt{

    // A simple scalar valued function interface, f : X -> R
    template <
        typename Real,
        template <typename> class XX
    >
    struct ScalarValuedFunction {
    private:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

    public:
        // <- f(x) 
        virtual Real operator () (const Vector& x) const = 0;

        // grad = grad f(x) 
        virtual void grad(const Vector& x,Vector& grad) const = 0;

        // H_dx = hess f(x) dx 
        virtual void hessvec(const Vector& x,const Vector& dx,Vector& H_dx)
            const = 0;

        // Allow a derived class to deallocate memory
        virtual ~ScalarValuedFunction() {}
    };
    
    template <
        typename Real,
        template <typename> class XX
    >
    struct ScalarValuedFunctionModifications {
    private:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

        // Disallow the copy constructor and the assignment operator
        ScalarValuedFunctionModifications(
            const ScalarValuedFunctionModifications&);
        ScalarValuedFunctionModifications& operator = (
            const ScalarValuedFunctionModifications&);

    public:
        // Use an empty default constructor
        ScalarValuedFunctionModifications() {}

        // Allow derived classes to deallocate memory
        virtual ~ScalarValuedFunctionModifications() {}

        // Merit function additions to the objective
        virtual Real merit(const Vector& x,const Real& f_x) const {
            return f_x;
        }

        // Stopping condition modification of the gradient
        virtual void grad_stop(
            const Vector& x,
            const Vector& grad,
            Vector& grad_stop
        ) const {
            X::copy(grad,grad_stop);
        }

        // Diagnostic modification of the gradient
        virtual void grad_diag(
            const Vector& x,
            const Vector& grad,
            Vector& grad_diag
        ) const {
            X::copy(grad,grad_diag);
        }

        // Modification of the gradient when finding a trial step
        virtual void grad_step(
            const Vector& x,
            const Vector& grad,
            Vector& grad_step
        ) const {
            X::copy(grad,grad_step);
        }

        // Modification of the gradient for a quasi-Newton method 
        virtual void grad_quasi(
            const Vector& x,
            const Vector& grad,
            Vector& grad_quasi
        ) const {
            X::copy(grad,grad_quasi);
        }

        // Modification of the gradient when solving for the equality multiplier
        virtual void grad_mult(
            const Vector& x,
            const Vector& grad,
            Vector& grad_mult
        ) const {
            X::copy(grad,grad_mult);
        }

        // Modification of the Hessian-vector product when finding a trial step
        virtual void hessvec_step(
            const Vector& x,
            const Vector& dx,
            const Vector& H_dx,
            Vector& Hdx_step 
        ) const {
            X::copy(H_dx,Hdx_step);
        }
    };


    // A simple vector valued function interface, f : X -> Y
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
        virtual void operator () (const X_Vector& x,Y_Vector& y) const = 0;

         // y=f'(x)dx 
         virtual void p(
             const X_Vector& x,
             const X_Vector& dx,
             Y_Vector& y
         ) const = 0;

         // z=f'(x)*dy
         virtual void ps(
             const X_Vector& x,
             const Y_Vector& dy,
             X_Vector& z
         ) const= 0;
         
         // z=(f''(x)dx)*dy
         virtual void pps(
             const X_Vector& x,
             const X_Vector& dx,
             const Y_Vector& dy,
             X_Vector& z
         ) const = 0;
         
         // Allow a derived class to deallocate memory
         virtual ~VectorValuedFunction() {}
    };

    // Defines how we output messages to the user
    struct Messaging {
        // Prints a message
        virtual void print(const std::string msg) const {
            std::cout << msg << std::endl;
        }

        // Prints an error
        virtual void error(const std::string msg) const {
            std::cerr << msg << std::endl;
            exit(EXIT_FAILURE);
        }

        // Allow a derived class to deallocate memory
        virtual ~Messaging() {}
    };

    // Which algorithm class do we use
    namespace AlgorithmClass{
        enum t{
            TrustRegion,            // Trust-Region algorithms
            LineSearch,             // Line-search algorithms
            UserDefined             // User provides the iterate 
        };

        // Converts the algorithm class to a string
        static std::string to_string(t algorithm_class){
            switch(algorithm_class){
            case TrustRegion:
                return "TrustRegion";
            case LineSearch:
                return "LineSearch";
            case UserDefined:
                return "UserDefined";
            default:
                throw;
            }
        }
        
        // Converts a string to an algorithm class 
        static t from_string(std::string algorithm_class){
            if(algorithm_class=="TrustRegion")
                return TrustRegion;
            else if(algorithm_class=="LineSearch")
                return LineSearch;
            else if(algorithm_class=="UserDefined")
                return UserDefined;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="TrustRegion" ||
                    name=="LineSearch" ||
                    name=="UserDefined"
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Reasons why we stop the algorithm
    namespace StoppingCondition{
        enum t{
            NotConverged,            // Algorithm did not converge
            RelativeGradientSmall,   // Relative gradient was sufficiently small
            RelativeStepSmall,       // Relative change in the step is small
            MaxItersExceeded,        // Maximum number of iterations exceeded
            InteriorPointInstability,// Instability in the interior point method
            UserDefined              // Some user defined stopping condition 
        };

        // Converts the stopping condition to a string 
        inline std::string to_string(t opt_stop){
            switch(opt_stop){
            case NotConverged:
                return "NotConverged";
            case RelativeGradientSmall:
                return "RelativeGradientSmall";
            case RelativeStepSmall:
                return "RelativeStepSmall";
            case MaxItersExceeded:
                return "MaxItersExceeded";
            case InteriorPointInstability:
                return "InteriorPointInstability";
            case UserDefined:
                return "UserDefined";
            default:
                throw;
            }
        }

        // Converts a string to a stopping condition
        inline t from_string(std::string opt_stop){
            if(opt_stop=="NotConverged")
                return NotConverged;
            else if(opt_stop=="RelativeGradientSmall")
                return RelativeGradientSmall;
            else if(opt_stop=="RelativeStepSmall")
                return RelativeStepSmall;
            else if(opt_stop=="MaxItersExceeded")
                return MaxItersExceeded;
            else if(opt_stop=="InteriorPointInstability")
                return InteriorPointInstability; 
            else if(opt_stop=="UserDefined")
                return UserDefined;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="NotConverged" ||
                    name=="RelativeGradientSmall" ||
                    name=="RelativeStepSmall" ||
                    name=="MaxItersExceeded" ||
                    name=="InteriorPointInstability" ||
                    name=="UserDefined"
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Various operators for both Hessian approximations and preconditioners
    namespace Operators{
        enum t{
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
        inline std::string to_string(t op){
            switch(op){
            case Identity:
                return "Identity";
            case ScaledIdentity:
                return "ScaledIdentity";
            case BFGS:
                return "BFGS";
            case InvBFGS:
                return "InvBFGS";
            case SR1:
                return "SR1";
            case InvSR1:
                return "InvSR1";
            case UserDefined:
                return "UserDefined";
            default:
                throw;
            }
        }
        
        // Converts a string to a operator 
        inline t from_string(std::string op){
            if(op=="Identity")
                return Identity; 
            else if(op=="ScaledIdentity")
                return ScaledIdentity; 
            else if(op=="BFGS")
                return BFGS; 
            else if(op=="InvBFGS")
                return InvBFGS; 
            else if(op=="SR1")
                return SR1; 
            else if(op=="InvSR1")
                return InvSR1; 
            else if(op=="UserDefined")
                return UserDefined; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="Identity" ||
                    name=="ScaledIdentity" ||
                    name=="BFGS" ||
                    name=="InvBFGS" ||
                    name=="SR1" ||
                    name=="InvSR1" ||
                    name=="UserDefined" 
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Different kinds of search directions 
    namespace LineSearchDirection{
        enum t{
            SteepestDescent,          // SteepestDescent 
            FletcherReeves,           // Fletcher-Reeves CG
            PolakRibiere,             // Polak-Ribiere CG
            HestenesStiefel,          // HestenesStiefel CG
            BFGS,                     // Limited-memory BFGS 
            NewtonCG                  // Newton-CG
        };
        
        // Converts the line-search direction to a string 
        inline std::string to_string(t dir){
            switch(dir){
            case SteepestDescent:
                return "SteepestDescent";
            case FletcherReeves:
                return "FletcherReeves";
            case PolakRibiere:
                return "PolakRibiere";
            case HestenesStiefel:
                return "HestenesStiefel";
            case BFGS:
                return "BFGS";
            case NewtonCG:
                return "NewtonCG";
            default:
                throw;
            }
        }
        
        // Converts a string to a line-search direction 
        inline t from_string(std::string dir){
            if(dir=="SteepestDescent")
                return SteepestDescent; 
            else if(dir=="FletcherReeves")
                return FletcherReeves; 
            else if(dir=="PolakRibiere")
                return PolakRibiere; 
            else if(dir=="HestenesStiefel")
                return HestenesStiefel; 
            else if(dir=="BFGS")
                return BFGS; 
            else if(dir=="NewtonCG")
                return NewtonCG; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="SteepestDescent" ||
                    name=="FletcherReeves" ||
                    name=="PolakRibiere" ||
                    name=="HestenesStiefel" ||
                    name=="BFGS" ||
                    name=="NewtonCG"
                )
                    return true;
                else
                    return false;
            }
        };
    };

    namespace LineSearchKind{
        enum t{
            Brents,           // Brent's minimization
            GoldenSection,    // Golden-section search 
            BackTracking,     // BackTracking search 
            TwoPointA,        // Barzilai and Borwein's method A
            TwoPointB         // Barzilai and Borwein's method B
        };
            
        // Converts the line-search kind to a string 
        inline std::string to_string(t kind){
            switch(kind){
            case Brents:
                return "Brents";
            case GoldenSection:
                return "GoldenSection";
            case BackTracking:
                return "BackTracking";
            case TwoPointA:
                return "TwoPointA";
            case TwoPointB:
                return "TwoPointB";
            default:
                throw;
            }
        }
        
        // Converts a string to a line-search kind 
        inline t from_string(std::string kind){
            if(kind=="Brents")
                return Brents; 
            else if(kind=="GoldenSection")
                return GoldenSection; 
            else if(kind=="BackTracking")
                return BackTracking; 
            else if(kind=="TwoPointA")
                return TwoPointA; 
            else if(kind=="TwoPointB")
                return TwoPointB; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="Brents" ||
                    name=="GoldenSection" ||
                    name=="BackTracking" ||
                    name=="TwoPointA" ||
                    name=="TwoPointB" 
                )
                    return true;
                else
                    return false;
            }
        };
    };
    
    namespace OptimizationLocation{
        enum t{
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
        inline std::string to_string(t loc){
            switch(loc){
            case BeginningOfOptimization:
                return "BeginningOfOptimization";
            case BeforeInitialFuncAndGrad:
                return "BeforeInitialFuncAndGrad";
            case AfterInitialFuncAndGrad:
                return "AfterInitialFuncAndGrad";
            case BeforeOptimizationLoop:
                return "BeforeOptimizationLoop";
            case BeforeSaveOld:
                return "BeforeSaveOld";
            case BeforeStep:
                return "BeforeStep";
            case BeforeGetStep:
                return "BeforeGetStep";
            case GetStep:
                return "GetStep";
            case AfterStepBeforeGradient:
                return "AfterStepBeforeGradient";
            case AfterGradient:
                return "AfterGradient";
            case BeforeQuasi:
                return "BeforeQuasi";
            case AfterQuasi:
                return "AfterQuasi";
            case EndOfOptimizationIteration:
                return "EndOfOptimizationIteration";
            case BeforeLineSearch:
                return "BeforeLineSearch";
            case AfterRejectedTrustRegion:
                return "AfterRejectedTrustRegion";
            case AfterRejectedLineSearch:
                return "AfterRejectedLineSearch";
            case BeforeActualVersusPredicted:
                return "BeforeActualVersusPredicted";
            case EndOfKrylovIteration:
                return "EndOfKrylovIteration";
            case EndOfOptimization:
                return "EndOfOptimization";
            default:
                throw;
            }
        }
        
        // Converts a string to a line-search kind 
        inline t from_string(std::string loc){
            if(loc=="BeginningOfOptimization")
                return BeginningOfOptimization;
            else if(loc=="BeforeInitialFuncAndGrad")
                return BeforeInitialFuncAndGrad;
            else if(loc=="AfterInitialFuncAndGrad")
                return AfterInitialFuncAndGrad;
            else if(loc=="BeforeOptimizationLoop")
                return BeforeOptimizationLoop;
            else if(loc=="BeforeSaveOld")
                return BeforeSaveOld; 
            else if(loc=="BeforeStep")
                return BeforeStep; 
            else if(loc=="BeforeGetStep")
                return BeforeGetStep; 
            else if(loc=="GetStep")
                return GetStep; 
            else if(loc=="AfterStepBeforeGradient")
                return AfterStepBeforeGradient; 
            else if(loc=="AfterGradient")
                return AfterGradient;
            else if(loc=="BeforeQuasi")
                return BeforeQuasi;
            else if(loc=="AfterQuasi")
                return AfterQuasi;
            else if(loc=="EndOfOptimizationIteration")
                return EndOfOptimizationIteration; 
            else if(loc=="BeforeLineSearch")
                return BeforeLineSearch;
            else if(loc=="AfterRejectedTrustRegion")
                return AfterRejectedTrustRegion;
            else if(loc=="AfterRejectedLineSearch")
                return AfterRejectedLineSearch;
            else if(loc=="BeforeActualVersusPredicted")
                return BeforeActualVersusPredicted;
            else if(loc=="EndOfKrylovIteration")
                return EndOfKrylovIteration;
            else if(loc=="EndOfOptimization")
                return EndOfOptimization;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="BeginningOfOptimization" ||
                    name=="BeforeInitialFuncAndGrad" ||
                    name=="AfterInitialFuncAndGrad" ||
                    name=="BeforeOptimizationLoop" ||
                    name=="BeforeSaveOld" || 
                    name=="BeforeStep" || 
                    name=="BeforeGetStep" || 
                    name=="GetStep" || 
                    name=="AfterStepBeforeGradient" ||
                    name=="AfterGradient" ||
                    name=="BeforeQuasi" ||
                    name=="AfterQuasi" ||
                    name=="EndOfOptimizationIteration" ||
                    name=="BeforeLineSearch" ||
                    name=="AfterRejectedTrustRegion" ||
                    name=="AfterRejectedLineSearch" ||
                    name=="BeforeActualVersusPredicted" ||
                    name=="EndOfKrylovIteration" ||
                    name=="EndOfOptimization"
                )
                    return true;
                else
                    return false;
            }
        };
    };
    
    // Different problem classes that we allow 
    namespace ProblemClass{
        enum t{
            Unconstrained,         // Unconstrained optimization 
            EqualityConstrained,   // Equality constrained optimization 
            InequalityConstrained, // Inequality constrained optimization 
            Constrained            // Fully constrained optimization 
        };

        // Converts the problem class to a string
        inline std::string to_string(t problem_class){
            switch(problem_class){
            case Unconstrained:
                return "Unconstrained";
            case EqualityConstrained:
                return "EqualityConstrained";
            case InequalityConstrained:
                return "InequalityConstrained";
            case Constrained:
                return "Constrained";
            default:
                    throw;
            }
        }

        // Converts a string to a problem class 
        inline t from_string(std::string problem_class){
            if(problem_class=="Unconstrained")
                return Unconstrained;
            else if(problem_class=="EqualityConstrained")
                return EqualityConstrained;
            else if(problem_class=="InequalityConstrained")
                return InequalityConstrained;
            else if(problem_class=="Constrained")
                return Constrained;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="Unconstrained" ||
                    name=="EqualityConstrained" ||
                    name=="InequalityConstrained" ||
                    name=="Constrained" 
                )
                    return true;
                else
                    return false;
            }
        };
    };
    
    // Different truncated Krylov solvers 
    namespace KrylovSolverTruncated{
        enum t{
            ConjugateDirection,         // Conjugate direction 
            MINRES                      // MINRES 
        };

        // Converts the problem class to a string
        inline std::string to_string(t truncated_krylov){
            switch(truncated_krylov){
            case ConjugateDirection:
                return "ConjugateDirection";
            case MINRES:
                return "MINRES";
            default:
                    throw;
            }
        }

        // Converts a string to a problem class 
        inline t from_string(std::string truncated_krylov){
            if(truncated_krylov=="ConjugateDirection")
                return ConjugateDirection;
            else if(truncated_krylov=="MINRES")
                return MINRES;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="ConjugateDirection" ||
                    name=="MINRES" 
                )
                    return true;
                else
                    return false;
            }
        };
    };
    
    
    // Different kinds of interior point methods
    namespace InteriorPointMethod{
        enum t{
            PrimalDual,          // Standard primal-dual interior point method 
            PrimalDualLinked,    // A primal dual IPM, but the primal and dual
                                 // variables are kept in lock step.
            LogBarrier           // Primal log-barrier method 
        };
        
        // Converts the interior point method to a string 
        inline std::string to_string(t ipm){
            switch(ipm){
            case PrimalDual:
                return "PrimalDual";
            case PrimalDualLinked:
                return "PrimalDualLinked";
            case LogBarrier:
                return "LogBarrier";
            default:
                throw;
            }
        }
        
        // Converts a string to an interior point method 
        inline t from_string(std::string ipm){
            if(ipm=="PrimalDual")
                return PrimalDual; 
            else if(ipm=="PrimalDualLinked")
                return PrimalDualLinked; 
            else if(ipm=="LogBarrier")
                return LogBarrier; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="PrimalDual" ||
                    name=="PrimalDualLinked" ||
                    name=="LogBarrier"
                )
                    return true;
                else
                    return false;
            }
        };
    };
    
    // Different schemes for adjusting the interior point centrality 
    namespace CentralityStrategy{
        enum t{
            Constant,           // We keep sigma fixed at each iteration.
            StairStep,          // If the relative improvement in the
                                // interior point parameter does not exceed
                                // that of the gradient, we do a constant
                                // reduction.  Otherwise, we hold sigma
                                // constant.
            PredictorCorrector  // On odd iterations, sigma=1, on even, sigma=0.
        };
        
        // Converts the centrality strategy to a string
        inline std::string to_string(t cstrat){
            switch(cstrat){
            case Constant:
                return "Constant";
            case StairStep:
                return "StairStep";
            case PredictorCorrector:
                return "PredictorCorrector";
            default:
                throw;
            }
        }
        
        // Converts a string to the cstrat
        inline t from_string(std::string cstrat){
            if(cstrat=="Constant")
                return Constant; 
            else if(cstrat=="StairStep")
                return StairStep; 
            else if(cstrat=="PredictorCorrector")
                return PredictorCorrector; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="Constant" ||
                    name=="StairStep" ||
                    name=="PredictorCorrector" 
                )
                    return true;
                else
                    return false;
            }
        };
    };


    // A collection of short routines that are only required locally
    namespace {
        // Returns the smallest positive non-Nan number between the two. 
        // If both are NaN, it will return NaN.
        template <typename Real>
        Real get_smallest(const Real x,const Real y) {
            return (x < y) || (y != y) ? x : y;
        }

        // Different nonlinear-CG directions 
        struct NonlinearCGDirections {
            enum t{
                HestenesStiefel,        
                PolakRibiere,           
                FletcherReeves          
            };
        };
            
        // Checks whether all the labels in the list labels are actually
        // labels stored in the is_label function.  If not, this function
        // throws an error.
        template <typename is_label>
        void checkLabels(
            const Messaging& msg,
            const std::list <std::string>& labels, 
            const std::string& kind
        ) {
                // Create a base message
                const std::string base
                    ="During serialization, found an invalid ";

                // Check the labels
                std::list <std::string>::const_iterator name = find_if(
                    labels.begin(), labels.end(),
                    std::not1(is_label()));

                if(name!=labels.end()) {
                    std::stringstream ss;
                    ss << base << kind << *name;
                    msg.error(ss.str());
                }
        }

        // Combines two strings in a funny sort of way.  Basically, given
        // a and b, if a is empty, this function returns b.  However, if
        // a is nonempty, it returns a.
        struct combineStrings : public std::binary_function
            <std::string,std::string,std::string>
        {
            std::string operator() (std::string a,std::string b) const {
                if(a=="") return b; else return a;
            }
        };
      
        // This checks the parameters and prints and error in there's a problem.
        template <typename checkParamVal>
        void checkParams(
            const Messaging& msg,
            const std::pair <std::list <std::string>,std::list <std::string> >&
                params 
        ) {
            std::string err=std::inner_product(
                params.first.begin(),
                params.first.end(),
                params.second.begin(),
                std::string(""),
                checkParamVal(),
                combineStrings()
            );
            if(err!="") msg.error(err);
        }

        // Converts a variety of basic datatypes to strings
        inline std::ostream& formatReal(std::ostream& out) {
            return out<<std::setprecision(2) << std::scientific << std::setw(10)
                << std::left;
        }
        inline std::ostream& formatInt(std::ostream& out) {
            return out << std::setw(10) << std::left;
        }
        inline std::ostream& formatString(std::ostream& out) {
            return out << std::setw(10) << std::left;
        }
        
        // Converts anything to a string.
        template <typename T>
        std::string atos(T x);
        template <>
        inline std::string atos <double> (double x){
            std::stringstream ss;
            ss << formatReal << x;
            return ss.str();
        }
        template <>
        inline std::string atos <Natural> (Natural x){
            std::stringstream ss;
            ss << formatInt << x;
            return ss.str();
        }
        template <>
        inline std::string atos <const char*> (const char* x){
            std::stringstream ss;
            ss << formatString << x;
            return ss.str();
        }
        template <>
        inline std::string atos <KrylovStop::t> (KrylovStop::t x){
            std::stringstream ss;
            // Converts the Krylov stopping condition to a shorter string 
            switch(x){
            case KrylovStop::NegativeCurvature:
                return atos <> ("NegCurv");
            case KrylovStop::RelativeErrorSmall:
                return atos <> ("RelErrSml");
            case KrylovStop::MaxItersExceeded:
                return atos <> ("IterExcd");
            case KrylovStop::TrustRegionViolated:
                return atos <> ("TrstReg");
            case KrylovStop::Instability:
                return atos <> ("Unstable");
            case KrylovStop::InvalidTrustRegionCenter:
                return atos <> ("InvldCnt");
            default:
                throw;
            }
        }

        // Blank separator for printing
        const std::string blankSeparator = ".         ";
    }

    // A collection of miscellaneous diagnostics that help determine errors.
    namespace Diagnostics {
        // Performs a 4-point finite difference directional derivative on
        // a scalar valued function f : X->R.  In other words, <- f'(x)dx.  We
        // accomplish this by doing a finite difference calculation on f.
        template <
            typename Real,
            template <typename> class XX
        >
        Real directionalDerivative(
            const ScalarValuedFunction<Real,XX>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const Real& epsilon
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx; X::init(x,x_op_dx);

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
            const ScalarValuedFunction<Real,XX>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const Real& epsilon,
            typename XX <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx; X::init(x,x_op_dx);

            // Create an element to store the gradient at this point 
            X_Vector fgrad_x_op_dx; X::init(x,fgrad_x_op_dx);

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
            const VectorValuedFunction<Real,XX,YY>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const Real& epsilon,
            typename YY <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx; X::init(x,x_op_dx);

            // Create an element for f(x+eps dx), etc.
            Y_Vector f_x_op_dx; Y::init(dd,f_x_op_dx);
            
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
            const VectorValuedFunction<Real,XX,YY>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const typename YY <Real>::Vector& dy,
            const Real& epsilon,
            typename XX <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx; X::init(x,x_op_dx);

            // Create an element for f'(x+eps dx)*dy, etc.
            X_Vector fps_xopdx_dy; X::init(dd,fps_xopdx_dy);

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
            const Messaging& msg,
            const ScalarValuedFunction<Real,XX>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Calculate the gradient at the point x
            X_Vector f_grad; X::init(x,f_grad);
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
            const Messaging& msg,
            const ScalarValuedFunction<Real,XX>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Create an element for the residual between the directional 
            // derivative computed Hessian-vector product and the true 
            // Hessian-vector product.
            X_Vector res; X::init(x,res);

            // Calculate hess f in the direction dx.  
            X_Vector hess_f_dx; X::init(x,hess_f_dx);
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
            const Messaging& msg,
            const ScalarValuedFunction<Real,XX>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const typename XX <Real>::Vector& dxx
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

            // Calculate hess f in the direction dx.  
            X_Vector H_x_dx; X::init(x,H_x_dx);
            f.hessvec(x,dx,H_x_dx);
            
            // Calculate hess f in the direction dxx.  
            X_Vector H_x_dxx; X::init(x,H_x_dxx);
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
            const Messaging& msg,
            const VectorValuedFunction<Real,XX,YY>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const typename YY <Real>::Vector& y
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            Y_Vector res; Y::init(y,res);

            // Calculate f'(x)dx 
            Y_Vector fp_x_dx; Y::init(y,fp_x_dx);
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
            const Messaging& msg,
            const VectorValuedFunction<Real,XX,YY>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const typename YY <Real>::Vector& dy
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Calculate f'(x)dx 
            Y_Vector fp_x_dx; Y::init(dy,fp_x_dx);
            f.p(x,dx,fp_x_dx);
            
            // Calculate f'(x)*dy 
            X_Vector fps_x_dy; X::init(dx,fps_x_dy);
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
            const Messaging& msg,
            const VectorValuedFunction<Real,XX,YY>& f,
            const typename XX <Real>::Vector& x,
            const typename XX <Real>::Vector& dx,
            const typename YY <Real>::Vector& dy
        ) {
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for the residual between the directional 
            // derivative and the true derivative.
            X_Vector res; X::init(x,res);

            // Calculate (f''(x)dx)*dy
            X_Vector fpps_x_dx_dy; X::init(dx,fpps_x_dx_dy);
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
    class StateManipulator {
    public:
        // Application
        virtual void operator () (
            const typename ProblemClass::Functions::t& fns,
            typename ProblemClass::State::t& state,
            OptimizationLocation::t loc
        ) const {};

        // Allow the derived class to deallocate memory
        virtual ~StateManipulator() {}
    };
   
    // A state manipulator that's been customized in order to print diagonistic
    // information
    template <typename ProblemClass>
    class DiagnosticManipulator : public StateManipulator <ProblemClass> {
    private:
        // A reference to an existing state manipulator 
        const StateManipulator <ProblemClass>& smanip;

        // A reference to the messsaging object
        const Messaging& msg;

    public:

        // Create a reference to an existing manipulator 
        explicit DiagnosticManipulator(
            const StateManipulator <ProblemClass>& smanip_,
            const Messaging& msg_
        ) : smanip(smanip_), msg(msg_) {}

        // Application
        void operator () (
            const typename ProblemClass::Functions::t& fns,
            typename ProblemClass::State::t& state,
            OptimizationLocation::t loc
        ) const {

            // Call the internal manipulator 
            smanip(fns,state,loc);

            // Create some shortcuts
            const Natural& msg_level=state.msg_level;

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
    // casting can file, so make sure to only use this when compatibility
    // can be guaranteed.
    template <typename Internal,typename External> 
    struct ConversionManipulator : public StateManipulator <External> {
    private:
        // A reference to the user-defined state manipulator
        const StateManipulator<Internal>& smanip;

    public:
        explicit ConversionManipulator(
            const StateManipulator <Internal>& smanip_
        ) : smanip(smanip_) {}

        // Application
        void operator () (
            const typename External::Functions::t& fns_,
            typename External::State::t& state_,
            OptimizationLocation::t loc
        ) const {
            const typename Internal::Functions::t& fns
                =dynamic_cast <const typename Internal::Functions::t&> (fns_);
            typename Internal::State::t& state 
                =dynamic_cast <typename Internal::State::t&> (state_);
            smanip(fns,state,loc);
        }
    };
       
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
    private:
        // This is a templated namespace.  Do not allow construction.
        Unconstrained();

    public:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        typedef std::pair < std::list <std::string>,
                            std::list <Real> > Reals;
        typedef std::pair < std::list <std::string>,
                            std::list <Natural> > Nats;
        typedef std::pair < std::list <std::string>,
                            std::list <std::string> > Params; 
        typedef std::pair < std::list <std::string>,
                            std::list <X_Vector> > X_Vectors;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            State();

        public:
            // The actual internal state of the optimization
            struct t {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                t& operator = (const t&);
                t(const t&);

            public:
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
                std::list <X_Vector> x; 
                
                // Gradient, possibly of the objective, possibly of the
                // Lagrangian.  It depends on the context.
                std::list <X_Vector> grad;
                
                // Trial step 
                std::list <X_Vector> dx;
                
                // Old optimization variable 
                std::list <X_Vector> x_old; 
                
                // Old gradient 
                std::list <X_Vector> grad_old;
                
                // Old trial step 
                std::list <X_Vector> dx_old;

                // Contains the prior iteration information for the
                // quasi-Newton operators
                std::list <X_Vector> oldY;
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

                // Trust-region parameter for checking whether a step has been
                // accepted
                Real eta2;

                // Actual reduction 
                Real ared;

                // Predicted reduction
                Real pred;

                // Number of rejected trust-region steps
                Natural rejected_trustregion;

                // ------------- LINE-SEARCH ------------- 

                // Line-search step length
                Real alpha;

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
                t() {
                    Unconstrained <Real,XX>::State::init_params(*this);
                }
                explicit t(const X_Vector& x) {
                    Unconstrained <Real,XX>::State::init_params(*this);
                    Unconstrained <Real,XX>::State::init_vectors(*this,x);
                }
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // This sets all of the parameters possible that don't require
            // special memory allocation such as variables.
            static void init_params_(t& state){
                state.eps_grad=Real(1e-8);
                state.eps_dx=Real(1e-8);
                state.stored_history=0;
                state.history_reset=5;
                state.iter=1;
                state.iter_max=10;
                state.opt_stop=StoppingCondition::NotConverged;
                state.krylov_iter=0;
                state.krylov_iter_max=10;
                state.krylov_iter_total=0;
                state.krylov_orthog_max=1;
                state.krylov_stop=KrylovStop::RelativeErrorSmall;
                state.krylov_rel_err=Real(0.);
                state.eps_krylov=Real(1e-2);
                state.krylov_solver
                    =KrylovSolverTruncated::ConjugateDirection;
                state.algorithm_class=AlgorithmClass::TrustRegion;
                state.PH_type=Operators::Identity;
                state.H_type=Operators::UserDefined;
                state.norm_gradtyp=std::numeric_limits<Real>::quiet_NaN();
                state.norm_dxtyp=std::numeric_limits<Real>::quiet_NaN();
                state.f_x=std::numeric_limits<Real>::quiet_NaN();
                state.f_xpdx=std::numeric_limits<Real>::quiet_NaN();
                state.msg_level=1;
                state.delta=Real(1.);
                state.eta1=Real(.1);
                state.eta2=Real(.9);
                state.ared=std::numeric_limits<Real>::quiet_NaN();
                state.pred=std::numeric_limits<Real>::quiet_NaN();
                state.rejected_trustregion=0;
                state.alpha=Real(1.);
                state.linesearch_iter=0;
                state.linesearch_iter_max=5;
                state.linesearch_iter_total=0;
                state.eps_ls=Real(1e-2);
                state.dir=LineSearchDirection::SteepestDescent;
                state.kind=LineSearchKind::GoldenSection;
            }
            static void init_params(t& state){
                Unconstrained <Real,XX>::State::init_params_(state);
            }

            // This initializes all the variables required for unconstrained
            // optimization.  
            static void init_vectors_(t& state,const X_Vector& x) {
                state.x.clear();
                    state.x.push_back(X_Vector());
                    X::init(x,state.x.front());
                    X::copy(x,state.x.front());
                state.grad.clear();
                    state.grad.push_back(X_Vector());
                    X::init(x,state.grad.front());
                state.dx.clear();
                    state.dx.push_back(X_Vector());
                    X::init(x,state.dx.front()); 
                state.x_old.clear();
                    state.x_old.push_back(X_Vector());
                    X::init(x,state.x_old.front());
                state.grad_old.clear();
                    state.grad_old.push_back(X_Vector());
                    X::init(x,state.grad_old.front()); 
                state.dx_old.clear();
                    state.dx_old.push_back(X_Vector());
                    X::init(x,state.dx_old.front()); 
            }
            static void init_vectors(t& state,const X_Vector& x) {
                Unconstrained <Real,XX>::State::init_vectors_(state,x);
            }

            // Initialize everything
            static void init(t& state,const X_Vector& x) {
                init_params(state);
                init_vectors(state,x);
            }

            // Check that we have a valid set of parameters.  
            static void check_(const Messaging& msg,const t& state) {
                   
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

                // Check that the trust-region radius is positive
                else if(state.delta<=Real(0.))
                    ss << "The trust-region radius must be positive: delta = "
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

                // Check that the line-search step length is positive 
                else if(state.alpha <= Real(0.)) 
                    ss << "The line-search step length must be positive: "
                        "alpha = " << state.alpha;
                
                // Check that the stopping tolerance for the line-search
                // methods is positive
                else if(state.eps_ls <= Real(0.)) 
                    ss << "The tolerance for the line-search stopping "
                        "condition must be positive: eps_ls = " << state.eps_ls;

                // If there's an error, print it
                if(ss.str()!="") msg.error(ss.str());
            }
            static void check(const Messaging& msg,const t& state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
            }
        };

        // Utilities for restarting the optimization
        struct Restart {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Restart();

        public:
            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
                        name == "alpha" || 
                        name == "eps_ls"
                    ) 
                        return true;
                    else
                        return false;
                }
            };

            // Checks whether we have a valid natural number label.
            struct is_nat : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_param : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                const Messaging& msg,
                const Reals& reals,
                const Nats& nats,
                const Params& params,
                const X_Vectors& xs
            ) {
                peopt::checkLabels <is_real> (msg,reals.first," real name: ");
                peopt::checkLabels <is_nat> (msg,nats.first," natural name: ");
                peopt::checkLabels <is_param>
                    (msg,params.first," paramater name: ");
                peopt::checkLabels <is_x> (msg,xs.first," variable name: ");
            }

            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string,std::string,std::string>
            {
                std::string operator() (
                    std::string label,
                    std::string val
                ) {
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
                typename State::t& state, 
                X_Vectors& xs
            ) {
                // Move the memory of all variables into the list 
                xs.first.push_back("x");
                xs.second.splice(xs.second.end(),state.x);
                xs.first.push_back("grad");
                xs.second.splice(xs.second.end(),state.grad);
                xs.first.push_back("dx");
                xs.second.splice(xs.second.end(),state.dx);
                xs.first.push_back("x_old");
                xs.second.splice(xs.second.end(),state.x_old);
                xs.first.push_back("grad_old");
                xs.second.splice(xs.second.end(),state.grad_old);
                xs.first.push_back("dx_old");
                xs.second.splice(xs.second.end(),state.dx_old);

                // Write out the quasi-Newton information with sequential names
                {Natural i=1;
                for(typename std::list<X_Vector>::iterator y=state.oldY.begin();
                    y!=state.oldY.end();
                    y=state.oldY.begin()
                ){
                    std::stringstream ss;
                    ss << "oldY_" << i;
                    xs.first.push_back(ss.str());
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
                    xs.first.push_back(ss.str());
                    xs.second.splice(xs.second.end(),state.oldS,s);
                }}
            }
            
            // Copy out all non-variables.  This includes reals, naturals,
            // and parameters
            static void stateToScalars(
                typename State::t& state, 
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                
                // Copy in all the real numbers 
                reals.first.push_back("eps_grad");
                reals.second.push_back(state.eps_grad);
                reals.first.push_back("eps_dx");
                reals.second.push_back(state.eps_dx);
                reals.first.push_back("krylov_rel_err");
                reals.second.push_back(state.krylov_rel_err);
                reals.first.push_back("eps_krylov");
                reals.second.push_back(state.eps_krylov);
                reals.first.push_back("norm_gradtyp");
                reals.second.push_back(state.norm_gradtyp);
                reals.first.push_back("norm_dxtyp");
                reals.second.push_back(state.norm_dxtyp);
                reals.first.push_back("f_x");
                reals.second.push_back(state.f_x);
                reals.first.push_back("f_xpdx");
                reals.second.push_back(state.f_xpdx);
                reals.first.push_back("delta");
                reals.second.push_back(state.delta);
                reals.first.push_back("eta1");
                reals.second.push_back(state.eta1);
                reals.first.push_back("eta2");
                reals.second.push_back(state.eta2);
                reals.first.push_back("ared");
                reals.second.push_back(state.ared);
                reals.first.push_back("pred");
                reals.second.push_back(state.pred);
                reals.first.push_back("alpha");
                reals.second.push_back(state.alpha);
                reals.first.push_back("eps_ls");
                reals.second.push_back(state.eps_ls);

                // Copy in all the natural numbers
                nats.first.push_back("stored_history");
                nats.second.push_back(state.stored_history);
                nats.first.push_back("history_reset");
                nats.second.push_back(state.history_reset);
                nats.first.push_back("iter");
                nats.second.push_back(state.iter);
                nats.first.push_back("iter_max");
                nats.second.push_back(state.iter_max);
                nats.first.push_back("krylov_iter");
                nats.second.push_back(state.krylov_iter);
                nats.first.push_back("krylov_iter_max");
                nats.second.push_back(state.krylov_iter_max);
                nats.first.push_back("krylov_iter_total");
                nats.second.push_back(state.krylov_iter_total);
                nats.first.push_back("krylov_orthog_max");
                nats.second.push_back(state.krylov_orthog_max);
                nats.first.push_back("msg_level");
                nats.second.push_back(state.msg_level);
                nats.first.push_back("rejected_trustregion");
                nats.second.push_back(state.rejected_trustregion);
                nats.first.push_back("linesearch_iter");
                nats.second.push_back(state.linesearch_iter);
                nats.first.push_back("linesearch_iter_max");
                nats.second.push_back(state.linesearch_iter_max);
                nats.first.push_back("linesearch_iter_total");
                nats.second.push_back(state.linesearch_iter_total);

                // Copy in all the parameters
                params.first.push_back("krylov_solver");
                params.second.push_back(
                    KrylovSolverTruncated::to_string(
                        state.krylov_solver));
                params.first.push_back("algorithm_class");
                params.second.push_back(
                    AlgorithmClass::to_string(state.algorithm_class));
                params.first.push_back("opt_stop");
                params.second.push_back(
                    StoppingCondition::to_string(state.opt_stop));
                params.first.push_back("krylov_stop");
                params.second.push_back(
                    KrylovStop::to_string(state.krylov_stop));
                params.first.push_back("H_type");
                params.second.push_back(
                    Operators::to_string(state.H_type));
                params.first.push_back("PH_type");
                params.second.push_back(
                    Operators::to_string(state.PH_type));
                params.first.push_back("dir");
                params.second.push_back(
                    LineSearchDirection::to_string(state.dir));
                params.first.push_back("kind");
                params.second.push_back(
                    LineSearchKind::to_string(state.kind));
            }

            // Copy in all variables.  This assumes that the quasi-Newton
            // information is being read in order.
            static void vectorsToState(
                typename State::t& state,
                X_Vectors& xs
            ) {
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
                    if(*name0=="x")
                        state.x.splice(state.x.end(),xs.second,x0);
                    else if(*name0=="grad")
                        state.grad.splice(state.grad.end(),xs.second,x0);
                    else if(*name0=="dx")
                        state.dx.splice(state.dx.end(),xs.second,x0);
                    else if(*name0=="x_old")
                        state.x_old.splice(state.x_old.end(),xs.second,x0);
                    else if(*name0=="grad_old")
                        state.grad_old.splice(
                            state.grad_old.end(),xs.second,x0);
                    else if(*name0=="dx_old")
                        state.dx_old.splice(state.dx_old.end(),xs.second,x0);
                    else if(name0->substr(0,5)=="oldY_")
                        state.oldY.splice(state.oldY.end(),xs.second,x0);
                    else if(name0->substr(0,5)=="oldS_")
                        state.oldS.splice(state.oldS.end(),xs.second,x0);

                    // Remove the string corresponding to the element just
                    // spliced if splicing occured.
                    if(xs.first.size() != xs.second.size())
                        xs.first.erase(name0);
                }
            }

            // Copy in all non-variables.  This includes reals, naturals,
            // and parameters
            static void scalarsToState(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
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
                    else if(*name=="alpha") state.alpha=*real;
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
                typename State::t& state,
                X_Vectors& xs,
                Reals& reals,
                Nats& nats,
                Params& params
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
                const Messaging& msg,
                typename State::t& state,
                X_Vectors& xs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs);

                // Check the strings used to represent parameters
                checkParams <checkParamVal> (msg,params);

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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Functions();

        public:
            // Actual storage of the functions required
            struct t{
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of auto_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                t& operator = (const t&);
                t(const t&);

            public:
                // Objective function
                std::auto_ptr <ScalarValuedFunction <Real,XX> > f;

                // Objective function modifications
                std::auto_ptr <ScalarValuedFunctionModifications <Real,XX> >
                    f_mod;

                // Preconditioner for the Hessian of the objective
                std::auto_ptr <Operator <Real,XX,XX> > PH;

                // Symmetric, positive definite operator that changes the
                // scaling of the trust-region.
                std::auto_ptr <Operator <Real,XX,XX> > TRS;

                // Initialize all of the pointers to null
                t() : f(NULL), PH(NULL), TRS(NULL) {}
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // The identity operator 
            struct Identity : public Operator <Real,XX,XX> {
                void operator () (const X_Vector& dx,X_Vector& result) const{
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
                const ScalarValuedFunctionModifications <Real,XX>& f_mod;

                // Current iterate
                const X_Vector& x;

                // Gradient of the objective 
                const X_Vector& grad;

                // Maximum size of the trust-region radius
                const Real& delta;
            public:
                ScaledIdentity(
                    const typename Functions::t& fns,
                    const typename State::t& state
                ) : f_mod(*(fns.f_mod)),
                    x(state.x.back()),
                    grad(state.grad.back()),
                    delta(state.delta)
                {};

                void operator () (const X_Vector& dx,X_Vector& result) const{
                    // Determine the norm of the gradient
                    X_Vector grad_step;
                        X::init(grad,grad_step);
                        f_mod.grad_step(x,grad,grad_step);
                    Real norm_grad=sqrt(X::innr(grad_step,grad_step));

                    // Copy in the direction and scale it
                    X::copy(dx,result);
                    X::scal(norm_grad/(Real(2.)*delta),result);
                }
            };

            // The BFGS Hessian approximation.  
            /* Note, the formula we normally see for BFGS denotes the inverse
                Hessian approximation.  This is not the inverse, but the true
                Hessian approximation. */ 
            class BFGS : public Operator <Real,XX,XX> {
            private:
                // Stored quasi-Newton information
                const std::list<X_Vector>& oldY;
                const std::list<X_Vector>& oldS;

                // Messaging device in case the quasi-Newton information is bad
                const Messaging& msg;
            public:
                BFGS(
                    const Messaging& msg_,
                    const typename State::t& state
                ) : oldY(state.oldY), oldS(state.oldS), msg(msg_) {};

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
                void operator () (const X_Vector& p, X_Vector& result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the BFGS Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work(oldY.size(),p);
                    for(typename std::list <X_Vector>::iterator w=work.begin();
                        w!=work.end();
                        w++
                    ) X::init(p,*w);

                    // If we have no vectors in our history, we return the
                    // direction
                    X::copy(p,result);
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
                        const X_Vector& si=*si_iter;
                        const X_Vector& yi=*yi_iter;
                        const X_Vector& Bisi=*Bisi_iter;

                        // Determine <Bisi,si>
                        Real inner_Bisi_si=X::innr(Bisi,si);

                        // Determine <yi,si>
                        Real inner_yi_si=X::innr(yi,si);

                        // Determine <si,Bip>
                        Real inner_si_Bip=X::innr(si,result);

                        // Determine <yi,p>
                        Real inner_yi_p=X::innr(yi,p);

                        // Determine -<si,Bip>/<Bisi,si> Bisi + Bip.  Store in
                        // Bip.  This will become B_{i+1}p.
                        X::axpy(-inner_si_Bip/inner_Bisi_si,Bisi,result);

                        // Determine <yi,p>/<yi,si> yi + w where we calculated w
                        // in the line above.  This completes the calculation of
                        // B_{i+1}p
                        X::axpy(inner_yi_p/inner_yi_si,yi,result);

                        // Check whether or not we've calculated B_{i+1}p for
                        // the last time
                        if(Bisi_iter==work.begin()) break;

                        // Begin the calculation of B_{i+1}sj
                        while(si_iter!=sj_iter){
                            // Add some additional references to the iterators 
                            const X_Vector& sj=*sj_iter;
                            X_Vector& Bisj=*Bisj_iter;

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
            /* The oldY and oldS lists have the same structure as the BFGS
            preconditioner. */
            class SR1 : public Operator <Real,XX,XX> {
            private:

                // Stored quasi-Newton information
                const std::list<X_Vector>& oldY;
                const std::list<X_Vector>& oldS;

                // Messaging device in case the quasi-Newton information is bad
                const Messaging& msg;
            public:
                SR1(
                    const Messaging& msg_,
                    const typename State::t& state
                ) : oldY(state.oldY), oldS(state.oldS), msg(msg_) {};
                
                // Operator interface
                void operator () (const X_Vector& p,X_Vector& result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the SR1 Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <X_Vector> work(oldY.size(),p);
                    for(typename std::list <X_Vector>::iterator w=work.begin();
                        w!=work.end();
                        w++
                    ) X::init(p,*w);

                    // If we have no vectors in our history, we return the 
                    // direction
                    X::copy(p,result);
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
                        const X_Vector& si=*si_iter;
                        const X_Vector& yi=*yi_iter;
                        const X_Vector& Bisi=*Bisi_iter;

                        // Determine <yi,p>
                        Real inner_yi_p=X::innr(yi,p);

                        // Determine <Bisi,p>
                        Real inner_Bisi_p=X::innr(Bisi,p);

                        // Determine <yi,si>
                        Real inner_yi_si=X::innr(yi,si);

                        // Determine <Bisi,si>
                        Real inner_Bisi_si=X::innr(Bisi,si);

                        // Determine (<yi,p>-<Bisi,p>) / (<y_i,s_i>-<Bisi,si>).
                        // Store in alpha
                        Real alpha=(inner_yi_p-inner_Bisi_p)
                            / (inner_yi_si-inner_Bisi_si);

                        // Determine alpha y_i + Bip.  Store in result (which
                        // accumulate Bip).
                        X::axpy(alpha,yi,result);

                        // Then, add -alpha*Bisi to this result
                        X::axpy(-alpha,Bisi,result);

                        // Check whether or not we've calculated B_{i+1}p for 
                        // the last time
                        if(Bisi_iter==work.begin()) break;

                        // Begin the calculation of B_{i+1}sj
                        while(si_iter!=sj_iter){
                            // Add some additional references to the iterators 
                            const X_Vector& sj=*sj_iter;
                            X_Vector& Bisj=*Bisj_iter;

                            // Determine <yi,sj>
                            Real inner_yi_sj=X::innr(yi,sj);

                            // Determine <Bisi,sj>
                            Real inner_Bisi_sj=X::innr(Bisi,sj);

                            // Determine (<yi,p>-<Bisi,p>)/(<y_i,s_i>-<Bisi,si>)
                            // Store in beta 
                            Real beta= (inner_yi_sj-inner_Bisi_sj) /
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
            /* The oldY list has the following structure
                oldY[0] = y_k = grad J(u_k) - grad J(u_{k-1})
                oldY[1] = y_{k-1} = grad J(u_{k-1}) - grad J(u_{k-2})
                The oldS list has the following structure
                oldS[0] = s_k = u_k - u_k{-1}
                oldS[1] = s_{k-1} = u_{k-1} - u_k{k-2} */
            class InvBFGS : public Operator <Real,XX,XX> {
            private:
                // Stored quasi-Newton information
                const std::list<X_Vector>& oldY;
                const std::list<X_Vector>& oldS;

                // Messaging device in case the quasi-Newton information is bad
                const Messaging& msg;
            public:
                InvBFGS(
                    const Messaging& msg_,
                    const typename State::t& state
                ) : oldY(state.oldY), oldS(state.oldS), msg(msg_) {};
                
                // Operator interface
                void operator () (const X_Vector& p,X_Vector& result) const{

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

                    // Before we begin computing, copy p to our result 
                    X::copy(p,result);

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
                        const X_Vector& y_k=*(y_iter++);
                        const X_Vector& s_k=*(s_iter++);
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
                        const X_Vector& s_k=*(--s_iter);
                        const X_Vector& y_k=*(--y_iter);

                        // beta=rho_i <y_i,result>
                        Real beta= rho[--i] * X::innr(y_k,result);

                        // result=  (alpha_i-beta) s_i + result
                        X::axpy(alpha[i]-beta,s_k,result);
                    }
                }
            };
            
            // The inverse SR1 operator.  
            /* In this definition, we take a shortcut and simply use the SR1
                Hessian approximation where we swap Y and S.  The oldY and oldS
                lists have the same structure as the BFGS operator. */
            class InvSR1 : public Operator <Real,XX,XX> {
            private:
                // Store the SR1 operator
                SR1 sr1;
            public:
                InvSR1(
                    const Messaging& msg,
                    const typename State::t& state
                ) : sr1(msg,state) {};
                void operator () (const X_Vector& p,X_Vector& result) const{
                    sr1(p,result);
                }
            };

            // A scalar valued function that overrides the Hessian if need be. 
            struct HessianAdjustedFunction
                : public peopt::ScalarValuedFunction <Real,XX>
            {
            private:

                // Hessian approximation
                std::auto_ptr <Operator <Real,XX,XX> > H;

                // Underlying function
                std::auto_ptr <peopt::ScalarValuedFunction <Real,XX> > f;

                // This forces derived classes to call the constructor that
                // depends on the state
                HessianAdjustedFunction() {}

            public:
                // The constructor determines whether we really need to build
                // a Hessian-vector product or if we use an internal
                // approximation
                HessianAdjustedFunction(
                    const Messaging& msg,
                    const typename State::t& state,
                    typename Functions::t& fns
                ) : f(fns.f) {
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
                 Real operator () (const X_Vector& x) const {
                    return (*f)(x);
                 }

                 // grad = grad f(x) 
                 void grad(const X_Vector& x,X_Vector& grad) const {
                        f->grad(x,grad);
                 }

                 // H_dx = hess f(x) dx 
                 // This actually computes the Hessian-vector product.  In 
                 // essence, we may want to use a Hessian approximation 
                 // provided by the optimization routines.  The following
                 // routine selects whether or not we use the hessvec 
                 // provided by the user.
                 virtual void hessvec(
                     const X_Vector& x,
                     const X_Vector& dx,
                     X_Vector& H_dx 
                 ) const {
                     if(H.get()!=NULL) 
                        (*H)(dx,H_dx);
                     else
                        f->hessvec(x,dx,H_dx);
                 }
            };

            // Check that all the functions are defined
            static void check(const Messaging& msg,const t& fns) {
                // Check that objective function exists 
                if(fns.f.get()==NULL)
                    msg.error("Missing an objective function definition.");
                
                // Check that objective function modifications exists 
                if(fns.f_mod.get()==NULL)
                    msg.error("Missing an objective function modification "
                        "definition.");
                
                // Check that a preconditioner exists 
                if(fns.PH.get()==NULL)
                    msg.error("Missing a preconditioner definition.");
            }

            // Initialize any missing functions for just unconstrained
            // optimization.
            static void init_(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
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
                        if(fns.PH.get()==NULL)
                            msg.error("An externally defined preconditioner "
                                "must be provided explicitly.");
                        break;
                    default:
                        msg.error("Not a valid Hessian approximation.");
                        break;
                }

                // If a trust-region operator has not been provided, use the
                // identity.
                if(fns.TRS.get()==NULL)
                    fns.TRS.reset(new Identity());

                // Check that all functions are defined (namely, the 
                // objective).
                check(msg,fns);

                // Modify the objective function if necessary
                fns.f.reset(new HessianAdjustedFunction(msg,state,fns));
            }

            // Initialize any missing functions 
            static void init(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
            ) {
                Unconstrained <Real,XX>::Functions::init_(msg,state,fns);
            }
        };

        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Printer();

        public:
            // Gets the header for the state information
            static void getStateHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) {

                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;

                // Basic information
                out.push_back(atos <> ("Iter"));
                out.push_back(atos <> ("f(x)"));
                out.push_back(atos <> ("merit(x)"));
                out.push_back(atos <> ("||grad||"));
                out.push_back(atos <> ("||dx||"));

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    out.push_back(atos <> ("KryIter"));
                    out.push_back(atos <> ("KryErr"));
                    out.push_back(atos <> ("KryWhy"));
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    out.push_back(atos <> ("LSIter"));
                }

                // In case we're using a trust-region method 
                if(algorithm_class==AlgorithmClass::TrustRegion) {
                    out.push_back(atos <> ("ared"));
                    out.push_back(atos <> ("pred"));
                    out.push_back(atos <> ("ared/pred"));
                }
            }

            // Combines all of the state headers
            static void getStateHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
            }

            // Gets the state information for output
            static void getState_(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {

                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.back();
                const X_Vector& dx=state.dx.back();
                const X_Vector& grad=state.grad.back();
                const Natural& iter=state.iter;
                const Real& f_x=state.f_x;
                const Natural& krylov_iter=state.krylov_iter;
                const Real& krylov_rel_err=state.krylov_rel_err;
                const KrylovStop::t& krylov_stop=state.krylov_stop;
                const Natural& linesearch_iter=state.linesearch_iter;
                const Real& ared=state.ared;
                const Real& pred=state.pred;
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;
                const Natural& rejected_trustregion=state.rejected_trustregion;

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
                X_Vector grad_diag;
                    X::init(grad,grad_diag);
                    f_mod.grad_diag(x,grad,grad_diag);
                Real norm_grad=sqrt(X::innr(grad_diag,grad_diag));


                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                if(!noiter)
                    out.push_back(atos <> (iter));
                else
                    out.push_back(atos <> ("*"));
                out.push_back(atos <> (f_x));
                out.push_back(atos <> (merit_x));
                out.push_back(atos <> (norm_grad));
                if(!opt_begin)
                    out.push_back(atos <> (norm_dx));
                else
                    out.push_back(blankSeparator);

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    if(!opt_begin) {
                        out.push_back(atos <> (krylov_iter));
                        out.push_back(atos <> (krylov_rel_err));
                        out.push_back(atos <> (krylov_stop));
                    } else 
                        for(Natural i=0;i<3;i++) out.push_back(blankSeparator);
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    if(!opt_begin)
                        out.push_back(atos <> (linesearch_iter));
                    else 
                        out.push_back(blankSeparator);
                }
                
                // In case we're using a trust-region method
                if(algorithm_class==AlgorithmClass::TrustRegion) {
                    if(!opt_begin) {
                        out.push_back(atos <> (ared));
                        out.push_back(atos <> (pred));
                        out.push_back(atos <> (ared/pred));
                    } else  
                        for(Natural i=0;i<3;i++) out.push_back(blankSeparator);
                }

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
            }

            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    out.push_back(atos <> ("KrySubItr"));
                    out.push_back(atos <> ("KryTotItr"));
                    out.push_back(atos <> ("KrySubErr"));
                }
            }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                // Create some shortcuts
                const Natural& krylov_iter=state.krylov_iter;
                const Natural& krylov_iter_total=state.krylov_iter_total;
                const Real& krylov_rel_err=state.krylov_rel_err;

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                out.push_back(atos <> (krylov_iter));
                out.push_back(atos <> (krylov_iter_total));
                out.push_back(atos <> (krylov_rel_err));
                
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
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Algorithms();

        public:
            // Checks a set of stopping conditions
            static StoppingCondition::t checkStop(
                const typename Functions::t& fns, 
                const typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.back();
                const X_Vector& grad=state.grad.back();
                const X_Vector& dx=state.dx.back();
                const Real& norm_gradtyp=state.norm_gradtyp;
                const Real& norm_dxtyp=state.norm_dxtyp;
                const Natural& iter=state.iter;
                const Natural& iter_max=state.iter_max;
                const Real& eps_grad=state.eps_grad;
                const Real& eps_dx=state.eps_dx;

                // Find both the norm of the gradient and the step
                X_Vector grad_stop;
                    X::init(grad,grad_stop);
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
                const ScalarValuedFunction <Real,XX>& f;

                // Objective modifications
                const ScalarValuedFunctionModifications <Real,XX>& f_mod;

                // Store a reference to the base of the Hessian-vector product
                const X_Vector& x;

                // Allocate memory for the Hessian modification
                mutable X_Vector H_dx;

            public:
                // Take in the objective and the base point during construction 
                HessianOperator(
                    const ScalarValuedFunction <Real,XX>& f_,
                    const ScalarValuedFunctionModifications <Real,XX>& f_mod_,
                    const X_Vector& x_)
                : f(f_), f_mod(f_mod_), x(x_) {
                    X::init(x,H_dx);
                }

                // Basic application
                void operator () (const X_Vector& dx,X_Vector &Hdx_step) const {
                    f.hessvec(x,dx,H_dx);
                    f_mod.hessvec_step(x,dx,H_dx,Hdx_step);
                }
            };
        
            // Checks whether we accept or reject a step
            static bool checkStep(
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& dx=state.dx.front();
                const X_Vector& grad=state.grad.front();
                const Real& eta1=state.eta1;
                const Real& eta2=state.eta2;
                const Real& f_x=state.f_x;
                const KrylovStop::t& krylov_stop=state.krylov_stop;
                Real& delta=state.delta;
                Real& ared=state.ared;
                Real& pred=state.pred;
                Real& f_xpdx=state.f_xpdx;
                
                // Allocate memory for temporaries that we need
                X_Vector x_p_dx; X::init(x,x_p_dx);

                // Determine x+dx 
                X::copy(dx,x_p_dx);
                X::axpy(Real(1.),x,x_p_dx);

                // Determine merit(x)
                Real merit_x = f_mod.merit(x,f_x);
                
                // Determine H(x)dx
                X_Vector H_dx;
                    X::init(x,H_dx);
                    f.hessvec(x,dx,H_dx);
                X_Vector Hdx_step;
                    X::init(x,Hdx_step);
                    f_mod.hessvec_step(x,dx,H_dx,Hdx_step);

                // Determine the gradient
                X_Vector grad_step;
                    X::init(x,grad_step);
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
                const Messaging& msg,
                const StateManipulator <Unconstrained <Real,XX> >& smanip,
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const Operator <Real,XX,XX>& PH=*(fns.PH);
                const Operator <Real,XX,XX>& TRS=*(fns.TRS);
                const Real& eps_dx=state.eps_dx;
                const Real& eps_krylov=state.eps_krylov;
                const Natural& krylov_iter_max=state.krylov_iter_max;
                const Natural& krylov_orthog_max=state.krylov_orthog_max;
                const Real& delta=state.delta;
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const Real& norm_dxtyp=state.norm_dxtyp;
                const KrylovSolverTruncated::t& krylov_solver
                    = state.krylov_solver;
                Natural& rejected_trustregion=state.rejected_trustregion;
                X_Vector& dx=state.dx.front();
                Natural& krylov_iter=state.krylov_iter;
                Natural& krylov_iter_total=state.krylov_iter_total;
                Real& krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                Natural& history_reset=state.history_reset;
                
                // Allocate some memory for the scaled trial step and the
                // trust-region center
                X_Vector x_tmp1; X::init(x,x_tmp1);
                X_Vector dx_cp; X::init(x,dx_cp);
                X_Vector grad_step; X::init(x,grad_step);
                X_Vector minus_grad; X::init(x,minus_grad);

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                X_Vector& dx=state.dx.front();

                // Determine the gradient for the step computation
                X_Vector grad_step; X::init(grad,grad_step);
                f_mod.grad_step(x,grad,grad_step);

                // We take the steepest descent direction
                X::copy(grad_step,dx);
                X::scal(Real(-1.),dx);
            }
    
            // Nonlinear Conjugate Gradient
            static void NonlinearCG(
                const NonlinearCGDirections::t dir,
                const typename Functions::t& fns,
                typename State::t& state
            ) {
            
                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& dx_old=state.dx_old.front();
                Natural& iter=state.iter;
                X_Vector& dx=state.dx.front();

                // Determine the gradient for the step computation
                X_Vector grad_step; X::init(grad,grad_step);
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

                    // Find -grad+beta*dx_old
                    X::copy(grad_step,dx);
                    X::scal(Real(-1.),dx);
                    X::axpy(beta,dx_old,dx);

                    // We don't ever check the strong-Wolfe conditions, so
                    // hard check that we have a descent direction
                    if(X::innr(dx,grad_step) > 0) X::scal(Real(-1.),dx);
                }
            }

            // Fletcher-Reeves CG search direction
            static Real FletcherReeves(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& grad_old=state.grad_old.front();

                // Determine the gradient for the step computation
                X_Vector grad_step;
                    X::init(grad,grad_step);
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector grad_old_step;
                    X::init(grad,grad_old_step);
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Return the momentum parameter
                return X::innr(grad_step,grad_step)
                    / X::innr(grad_old_step,grad_old_step);
            }
        
            // Polak-Ribiere CG search direction
            static Real PolakRibiere(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& grad_old=state.grad_old.front();

                // Determine the gradient for the step computation
                X_Vector grad_step; X::init(grad,grad_step);
                X_Vector grad_old_step; X::init(grad,grad_old_step);
                f_mod.grad_step(x,grad,grad_step);
                f_mod.grad_step(x,grad_old,grad_old_step);

                // Find grad-grad_old 
                X_Vector grad_m_gradold; X::init(grad,grad_m_gradold);
                X::copy(grad_step,grad_m_gradold);
                X::axpy(Real(-1.),grad_old_step,grad_m_gradold);
                    
                // Return the momentum parameter
                return X::innr(grad_step,grad_m_gradold)
                    / X::innr(grad_old_step,grad_old_step);
            }
            
            // Hestenes-Stiefel search direction
            static Real HestenesStiefel(
                const typename Functions::t& fns,
                typename State::t& state
            ) {

                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& grad_old=state.grad_old.front();
                const X_Vector& dx_old=state.dx_old.front();
                const Real& alpha=state.alpha;

                // Determine the gradient for the step computation
                X_Vector grad_step;
                    X::init(grad,grad_step);
                    f_mod.grad_step(x,grad,grad_step);
                X_Vector grad_old_step;
                    X::init(grad,grad_old_step);
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Find grad-grad_old 
                X_Vector grad_m_gradold; X::init(grad,grad_m_gradold);
                X::copy(grad_step,grad_m_gradold);
                X::axpy(Real(-1.),grad_old_step,grad_m_gradold);
                    
                // Return the momentum parameter.  Note, we scale things by
                // alpha here, which is not in the standard Hestenes-Stiefel
                // formula.  Internal to this code, we scale our search
                // directions by alpha at every iteration.  This allows us
                // to use the update x+dx for both trust-region and line-search
                // methods.  Now, the formulas for Hestenes-Stiefel assume
                // that the directions have not been scaled.  Therefore,
                // dx_old really needs to be (1/alpha)dx_old to get rid of
                // the scaling.  The alpha on top is just an algebraic
                // rearrangement.
                return alpha*X::innr(grad_step,grad_m_gradold)
                    / X::innr(dx_old,grad_m_gradold);
            }

            // BFGS search direction
            static void BFGS(
                const Messaging& msg,
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                
                // Create some shortcuts 
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                X_Vector& dx=state.dx.front();

                // Determine the gradient for the step computation
                X_Vector grad_step; X::init(grad,grad_step);
                f_mod.grad_step(x,grad,grad_step);

                // Create the inverse BFGS operator
                typename Functions::InvBFGS Hinv(msg,state); 

                // Apply the inverse BFGS operator to the gradient
                Hinv(grad_step,dx);

                // Negate the result
                X::scal(Real(-1.),dx);
            }

            // Compute a Golden-Section search between eps and 2*alpha where
            // alpha is the last line search parameter.
            static void goldenSection(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                const X_Vector& x=state.x.front();
                const Natural& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                X_Vector& dx=state.dx.front();
                Natural& iter_total=state.linesearch_iter_total;
                Natural& iter=state.linesearch_iter;
                Real& f_xpdx=state.f_xpdx;
                
                // Create one work element that holds x+mu dx or x+lambda dx 
                X_Vector x_p_dx; X::init(x,x_p_dx);

                // Find 1 over the golden ratio
                Real beta=Real(2./(1.+sqrt(5.)));

                // Find a bracket for the linesearch such that a < b
                Real a=Real(0.);
                Real b=Real(2.)*alpha;

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

                // Search for a fixed number of iterations 
                for(iter=0;iter<iter_max;iter++){

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
                // guess for the line search parameter
                alpha=merit_lambda < merit_mu ? lambda : mu;

                // Save the objective value at this step
                f_xpdx=merit_lambda < merit_mu ? f_lambda : f_mu;
            }

            // Find the line search parameter based on the 2-point approximation
            // from Barzilai and Borwein
            static void twoPoint(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& x_old=state.x_old.front();
                const X_Vector& grad_old=state.grad_old.front();
                const LineSearchKind::t& kind=state.kind;
                Real& alpha=state.alpha;
                X_Vector& dx=state.dx.front();
                Natural& iter_total=state.linesearch_iter_total;
                Natural& iter=state.linesearch_iter;
                Real& f_xpdx=state.f_xpdx;
                
                // Create elements for delta_x and delta_grad as well as one
                // work element for storing x+alpha dx 
                X_Vector delta_x; X::init(x,delta_x);
                X_Vector delta_grad; X::init(x,delta_grad);
                X_Vector x_p_dx; X::init(x,x_p_dx);

                // Find delta_x
                X::copy(x,delta_x);
                X::axpy(Real(-1.),x_old,delta_x);

                // Determine the gradient for the step computation
                X_Vector grad_step;
                    X::init(grad,grad_step);
                    f_mod.grad_step(x,grad,grad_step);
                
                X_Vector grad_old_step;
                    X::init(grad,grad_old_step);
                    f_mod.grad_step(x,grad_old,grad_old_step);

                // Find delta_grad
                X::copy(grad_step,delta_grad);
                X::axpy(Real(-1.),grad_old_step,delta_grad);

                // Find alpha
                if(kind==LineSearchKind::TwoPointA)
                    alpha=X::innr(delta_x,delta_grad) 
                        / X::innr(delta_grad,delta_grad);
                else if(kind==LineSearchKind::TwoPointB)
                    alpha=X::innr(delta_x,delta_x)/X::innr(delta_x,delta_grad);

                // Save the objective value at this step
                X::copy(x,x_p_dx);
                X::axpy(alpha,dx,x_p_dx);
                f_xpdx=f(x_p_dx);

                // Since we do one function evaluation, increase the linesearch
                // iteration by one
                iter=1; iter_total++;
            }
            
            // Compute a backtracking line-search. 
            static void backTracking(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod 
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const Natural& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                X_Vector& dx=state.dx.front();
                Natural& iter_total=state.linesearch_iter_total;
                Natural& iter=state.linesearch_iter;
                Real& f_xpdx=state.f_xpdx;
                
                // Create one work element for holding x+alpha s
                X_Vector x_p_dx; X::init(x,x_p_dx);

                // Store the best merit value and alpha that we used to find it.
                // Our initial guess will be at alpha*2.
                Real alpha_best=Real(2.)*alpha;
                X::copy(x,x_p_dx);
                X::axpy(alpha_best,dx,x_p_dx);
                Real f_best=f(x_p_dx);
                Real merit_best=f_mod.merit(x_p_dx,f_best);

                // Evaluate the merit iter_max times at a distance of
                // 2*alpha, alpha, alpha/2, ....  Then, pick the best one.
                // Note, we start iter at 1 since we've already done one
                // iteration above.
                Real alpha0=alpha;
                for(iter=1;iter<iter_max;iter++){
                    // Evaluate f(x+alpha*dx)
                    X::copy(x,x_p_dx);
                    X::axpy(alpha0,dx,x_p_dx);
                    f_xpdx=f(x_p_dx);
                    Real merit_xpdx=f_mod.merit(x_p_dx,f_xpdx);

                    // If this is better than our best guess so far, save it
                    if(merit_xpdx < merit_best){
                        f_best=f_xpdx;
                        merit_best=merit_xpdx;
                        alpha_best=alpha0;
                    }

                    // Reduce the size of alpha
                    alpha0 /= Real(2.);
                }

                // Save the best merit value and alpha found
                alpha=alpha_best;
                f_xpdx=f_best;

                // Indicate how many iterations we used to find this value
                iter_total+=iter;
            }

            // Finds a trial step using a line-search for globalization
            static void getStepLS(
                const Messaging& msg,
                const StateManipulator <Unconstrained <Real,XX> >& smanip,
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const Operator <Real,XX,XX>& PH=*(fns.PH);
                const Operator <Real,XX,XX>& TRS=*(fns.TRS);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const LineSearchDirection::t& dir=state.dir;
                const LineSearchKind::t& kind=state.kind;
                const Natural& iter=state.iter;
                const Natural& linesearch_iter_max=state.linesearch_iter_max;
                const Real& f_x=state.f_x;
                const Real& eps_dx=state.eps_dx;
                const Real& norm_dxtyp=state.norm_dxtyp;
                const Real& eps_krylov=state.eps_krylov;
                const Natural& krylov_iter_max=state.krylov_iter_max;
                const Natural& krylov_orthog_max=state.krylov_orthog_max;
                const KrylovSolverTruncated::t& krylov_solver
                    = state.krylov_solver;
                X_Vector& dx=state.dx.front();
                Real& f_xpdx=state.f_xpdx;
                Real& alpha=state.alpha;
                Real& krylov_rel_err=state.krylov_rel_err;
                Natural& krylov_iter=state.krylov_iter;
                Natural& krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                
                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeGetStep);

                // Create the trust-region center 
                X_Vector x_cntr; X::init(x,x_cntr);
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
                    X_Vector dx_cp; X::init(x,dx_cp);
                    X_Vector grad_step;
                        X::init(grad,grad_step);
                        f_mod.grad_step(x,grad,grad_step);
                    X_Vector minus_grad;
                        X::init(x,minus_grad);
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

                // Do a line-search in the specified direction
                X_Vector x_p_dx; X::init(x,x_p_dx);
                Real merit_x(std::numeric_limits <Real>::quiet_NaN());
                Real merit_xpdx(std::numeric_limits <Real>::quiet_NaN());
                switch(kind){
                case LineSearchKind::GoldenSection:
                    // Continue doing a line-search until we get a reduction
                    // in the merit value.
                    do {
                        // Conduct the golden section search
                        goldenSection(fns,state);

                        // Find x+alpha dx
                        X::copy(x,x_p_dx);
                        X::axpy(alpha,dx,x_p_dx);

                        // If we have no reduction in the merit, print
                        // some diagnostic information.
                        merit_x = f_mod.merit(x,f_x);
                        merit_xpdx = f_mod.merit(x_p_dx,f_xpdx);
                        if(merit_xpdx > merit_x || merit_xpdx!=merit_xpdx) {

                            // Determine the size of the step
                            Real norm_dx=alpha*sqrt(X::innr(dx,dx));

                            // Check if the step becomes so small that we're not
                            // making progress.  In this case, take a zero step 
                            // and allow the stopping conditions to exit
                            if(norm_dx < eps_dx*norm_dxtyp) {
                                alpha=0.;
                                break;
                            }

                            // Manipulate the state if required
                            smanip(fns,state,
                                OptimizationLocation::AfterRejectedLineSearch);

                            // We reduce alpha by a factor of four when we
                            // reject the step since the line-search always
                            // looks twice alpha out in the direction of the 
                            // search direction.  By reducing alpha by a factor
                            // of four we insure that the next line-search
                            // examines a unique set of points.
                            alpha /= Real(4.);
                        }

                    // If we don't decrease the merit , try again 
                    } while(merit_x < merit_xpdx || merit_xpdx!=merit_xpdx);
                    break;
                case LineSearchKind::BackTracking:
                    // Continue doing a line-search until we get a reduction
                    // in the merit value.
                    do {
                        // Conduct the backtracking search
                        backTracking(fns,state);

                        // Find x+alpha dx
                        X::copy(x,x_p_dx);
                        X::axpy(alpha,dx,x_p_dx);

                        // If we have no reduction in the merit, print
                        // some diagnostic information.
                        merit_x = f_mod.merit(x,f_x);
                        merit_xpdx = f_mod.merit(x_p_dx,f_xpdx);
                        if(merit_xpdx > merit_x || merit_xpdx!=merit_xpdx) {
                            // Determine the size of the step
                            Real norm_dx=alpha*sqrt(X::innr(dx,dx));

                            // Check if the step becomes so small that we're not
                            // making progress.  In this case, take a zero step 
                            // and allow the stopping conditions to exit
                            if(norm_dx < eps_dx*norm_dxtyp) {
                                alpha = Real(0.);
                                break;
                            }

                            // Manipulate the state if required
                            smanip(fns,state,
                                OptimizationLocation::AfterRejectedLineSearch);

                            // We set alpha to be four times less than the
                            // minimimum alpha we searched before.  We do this
                            // since the line-search always looks twice alpha
                            // out in the beginning of the search.  In addition,
                            // we use a loop here instead of the pow routine
                            // since the cmath pow routine requires an integer
                            // second argument and, depending on the
                            // architecture, our internally defined integer
                            // may be too large for the routine.  This is
                            // unlikely, but some versions of gcc complain.
                            for(Natural i=1;i<=linesearch_iter_max+1;i++)
                                alpha = alpha/Real(2.);
                        }

                    // If we don't decrease the merit, try again 
                    } while(merit_x < merit_xpdx || merit_xpdx!=merit_xpdx);
                    break;
                case LineSearchKind::TwoPointA:
                case LineSearchKind::TwoPointB:
                    if(iter>1)
                        twoPoint(fns,state);
                    else
                        goldenSection(fns,state);
                    break;
                case LineSearchKind::Brents:
                    msg.error(
                        "Brent's linesearch is not currently implemented.");
                    break;
                }
            
                // Scale the line-search direction by the line search parameter 
                X::scal(alpha,dx);
            }

            // Finds a new trial step
            static void getStep(
                const Messaging& msg,
                const StateManipulator <Unconstrained <Real,XX> >& smanip,
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;

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
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Exit immediately if we're not using a quasi-Newton method
                if(state.stored_history==0) return;

                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& x_old=state.x_old.front();
                const X_Vector& grad_old=state.grad_old.front();
                const Operators::t& PH_type=state.PH_type;
                const Operators::t& H_type=state.H_type;
                const LineSearchDirection::t& dir=state.dir;
                std::list <X_Vector>& oldY=state.oldY;
                std::list <X_Vector>& oldS=state.oldS;
               
                // Allocate some temp storage for y and s
                X_Vector s; X::init(x,s);
                X_Vector y; X::init(x,y);

                // Find s = x-x_old
                X::copy(x,s);
                X::axpy(Real(-1.),x_old,s);
                
                // Determine the gradient for the quasi-Newton computation 
                X_Vector grad_quasi;
                    X::init(grad,grad_quasi);
                    f_mod.grad_quasi(x,grad,grad_quasi);
                X_Vector grad_old_quasi;
                    X::init(grad_old,grad_old_quasi);
                    f_mod.grad_quasi(x,grad_old,grad_old_quasi);

                // Find y = grad - grad_old
                X::copy(grad_quasi,y);
                X::axpy(Real(-1.),grad_old_quasi,y);

                // If we're using BFGS, check that <y,s> > 0
                if((PH_type==Operators::InvBFGS || H_type==Operators::BFGS
                    || dir==LineSearchDirection::BFGS)
                    && X::innr(y,s) <= Real(0.))
                    return;

                // Insert these into the quasi-Newton storage
                oldS.push_front(s);
                oldY.push_front(y);

                // Determine if we need to free some memory
                if(oldS.size()>state.stored_history){
                    oldS.pop_back();
                    oldY.pop_back();
                }
            }

            // Solves an optimization problem
            static void getMin_(
                const Messaging& msg,
                const StateManipulator <Unconstrained <Real,XX> >& smanip,
                typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>&
                    f_mod=*(fns.f_mod);
                X_Vector& x=state.x.front();
                X_Vector& grad=state.grad.front();
                X_Vector& dx=state.dx.front();
                X_Vector& x_old=state.x_old.front();
                X_Vector& grad_old=state.grad_old.front();
                X_Vector& dx_old=state.dx_old.front();
                Real& f_x=state.f_x;
                Real& f_xpdx=state.f_xpdx;
                Real& norm_gradtyp=state.norm_gradtyp;
                Real& norm_dxtyp=state.norm_dxtyp;
                Natural& iter=state.iter;
                StoppingCondition::t& opt_stop=state.opt_stop;
                
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
                    X_Vector grad_stop;
                        X::init(grad,grad_stop);
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
                const Messaging& msg,
                typename Functions::t& fns,
                typename State::t& state
            ){
                // Create an empty state manipulator
                StateManipulator <Unconstrained <Real,XX> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                const Messaging& msg,
                const StateManipulator <Unconstrained <Real,XX> >& smanip,
                typename Functions::t& fns,
                typename State::t& state
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
    private:
        // This is a templated namespace.  Do not allow construction.
        EqualityConstrained();

    public:
        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector;
        
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

        // This defines a product space between X and Y
        template <typename Real_>
        struct XXxYY {
            typedef std::pair <X_Vector,Y_Vector> Vector;

            // Memory allocation and size setting
            static void init(const Vector& x, Vector& y) {
                X::init(x.first,y.first);
                Y::init(x.second,y.second);
            }

            // y <- x (Shallow.  No memory allocation.)
            static void copy(const Vector& x, Vector& y) {
                X::copy(x.first,y.first);
                Y::copy(x.second,y.second);
            }

            // x <- alpha * x
            static void scal(const Real_& alpha, Vector& x) {
                X::scal(alpha,x.first);
                Y::scal(alpha,x.second);
            }

            // x <- 0 
            static void zero(Vector& x) {
                X::zero(x.first);
                Y::zero(x.second);
            }

            // y <- alpha * x + y
            static void axpy(const Real_& alpha, const Vector& x, Vector& y) {
                X::axpy(alpha,x.first,y.first);
                Y::axpy(alpha,x.second,y.second);
            }

            // innr <- <x,y>
            static Real_ innr(const Vector& x,const Vector& y) {
                return X::innr(x.first,y.first) + Y::innr(x.second,y.second);
            }
        };
        typedef XXxYY <Real> XxY;
        typedef typename XxY::Vector XxY_Vector;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            State();

        public:
            // The actual internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                t& operator = (const t&);
                t(const t&);

            public:
                // The Lagrange multiplier (dual variable) for the equality
                // constraints
                std::list <Y_Vector> y;

                // The Lagrange multiplier step
                std::list <Y_Vector> dy;

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
                Real xi_lmh;    // Lagrange multiplier

                // Sets all the inexactness tolerances 
                void xi_all(const Real& xi) {
                    xi_qn = xi;
                    xi_pg = xi;
                    xi_proj = xi;
                    xi_tang = xi;
                    xi_lmh = xi;
                }

                // Absolute tolerance on the residual of the Lagrange
                // multiplier solve.
                Real xi_lmg;

                // Tolerance for how much error is acceptable after computing
                // the tangential step given the result from the tangential
                // subproblem;
                Real xi_4;

                // Residual term in the predicted reduction
                Real rpred;

                // Preconditioners for the augmented system
                Operators::t PSchur_left_type;
                Operators::t PSchur_right_type;

                // Maximum number of iterations used when solving the augmented
                // system
                Natural augsys_iter_max;

                // How often we restart the augmented system solve
                Natural augsys_rst_freq;
                
                // Equality constraint evaluated at x.  This is used in the
                // quasinormal step as well as in the computation of the
                // linear Taylor series at x in the direciton dx_n.
                std::list <Y_Vector> g_x;

                // A typical norm for norm_gx.  Generally, we just take
                // the value at the first iteration.
                Real norm_gxtyp;
                
                // Linear Taylor series at x in the direction dx_n.  This is
                // used both in the predicted reduction as well as the
                // residual predicted reduction. 
                std::list <Y_Vector> gpxdxn_p_gx;

                // Derivative of the constraint applied to the tangential step
                // this is used in the residual predicted reduction.
                std::list <Y_Vector> gpxdxt;

                // Norm of gpxdxn_p_gx.  This is used in the penalty parameter
                // computation and predicted reduction. 
                Real norm_gpxdxnpgx;

                // Normal step
                std::list <X_Vector> dx_n;
                
                // Cauchy point for normal step
                std::list <X_Vector> dx_ncp;

                // (Corrected) tangential step
                std::list <X_Vector> dx_t;

                // Tangential step prior to correction
                std::list <X_Vector> dx_t_uncorrected;
                
                // Cauchy point for tangential step prior to correction
                std::list <X_Vector> dx_tcp_uncorrected;
                
                // Hessian applied to the normal step.  This is required by
                // W_gradpHdxn as well as the predicted reduction.
                std::list <X_Vector> H_dxn;

                // Quantity grad f(x) + g'(x)*y + H dx_n.  This is required
                // in the tangential subproblem and the predicted reduction.
                std::list <X_Vector> W_gradpHdxn;
                
                // Hessian applied to the uncorrected tangential step.  THis
                // is needed in the predicted reduction.
                std::list <X_Vector> H_dxtuncorrected;
                
                // Initialization constructors
                t() {
                    EqualityConstrained <Real,XX,YY>::State::init_params(*this);
                }
                t(const X_Vector& x,const Y_Vector& y) {
                    EqualityConstrained <Real,XX,YY>::State::init_params(*this);
                    EqualityConstrained <Real,XX,YY>::State
                        ::init_vectors(*this,x,y);
                }
            };
            
            // This initializes all the parameters required for equality
            // constrained optimization.  
            static void init_params_(t& state) {
                state.zeta = Real(0.8);
                state.eta0 = Real(0.5);
                state.rho = Real(1.0);
                state.rho_old = state.rho;
                state.rho_bar = Real(1e-8);
                state.eps_constr = Real(1e-8);
                state.xi_all(Real(1e-4));
                state.xi_lmg = Real(1e4);
                state.xi_4 = Real(2.);
                state.rpred=std::numeric_limits<Real>::quiet_NaN();
                state.norm_gxtyp=std::numeric_limits<Real>::quiet_NaN();
                state.norm_gpxdxnpgx=std::numeric_limits<Real>::quiet_NaN();
                state.PSchur_left_type=Operators::Identity;
                state.PSchur_right_type=Operators::Identity;
                state.augsys_iter_max = 100;
                state.augsys_rst_freq = 0;
            }
            static void init_params(t& state) {
                Unconstrained <Real,XX>::State::init_params_(state); 
                EqualityConstrained <Real,XX,YY>::State::init_params_(state);
            }

            // This initializes all the variables required for equality
            // constrained optimization.  
            static void init_vectors_(
                t& state,
                const X_Vector& x,
                const Y_Vector& y
            ) {
                state.y.clear();
                    state.y.push_back(Y_Vector());
                    Y::init(y,state.y.front());
                    Y::copy(y,state.y.front());
                
                state.dy.clear();
                    state.dy.push_back(Y_Vector());
                    Y::init(y,state.dy.front());
                
                state.g_x.clear();
                    state.g_x.push_back(Y_Vector());
                    Y::init(y,state.g_x.front());
                
                state.gpxdxn_p_gx.clear();
                    state.gpxdxn_p_gx.push_back(Y_Vector());
                    Y::init(y,state.gpxdxn_p_gx.front());
                
                state.gpxdxt.clear();
                    state.gpxdxt.push_back(Y_Vector());
                    Y::init(y,state.gpxdxt.front());
                
                state.dx_n.clear();
                    state.dx_n.push_back(X_Vector());
                    X::init(x,state.dx_n.front());
                
                state.dx_ncp.clear();
                    state.dx_ncp.push_back(X_Vector());
                    X::init(x,state.dx_ncp.front());
                
                state.dx_t.clear();
                    state.dx_t.push_back(X_Vector());
                    X::init(x,state.dx_t.front());
                
                state.dx_t_uncorrected.clear();
                    state.dx_t_uncorrected.push_back(X_Vector());
                    X::init(x,state.dx_t_uncorrected.front());
                
                state.dx_tcp_uncorrected.clear();
                    state.dx_tcp_uncorrected.push_back(X_Vector());
                    X::init(x,state.dx_tcp_uncorrected.front());
                
                state.H_dxn.clear();
                    state.H_dxn.push_back(X_Vector());
                    X::init(x,state.H_dxn.front());
                
                state.W_gradpHdxn.clear();
                    state.W_gradpHdxn.push_back(X_Vector());
                    X::init(x,state.W_gradpHdxn.front());
                
                state.H_dxtuncorrected.clear();
                    state.H_dxtuncorrected.push_back(X_Vector());
                    X::init(x,state.H_dxtuncorrected.front());
            }
            static void init_vectors(
                t& state,
                const X_Vector& x,
                const Y_Vector& y
            ) {
                Unconstrained <Real,XX>::State::init_vectors_(state,x); 
                EqualityConstrained <Real,XX,YY>::State
                    ::init_vectors_(state,x,y);
            }
           
            // Initializes everything
            static void init(t& state, const X_Vector& x, const Y_Vector& y) {
                init_params(state);
                init_vectors(state,x,y);
            }

            // Check that we have a valid set of parameters.  
            static void check_(const Messaging& msg,const t& state) {
                   
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
            static void check(const Messaging& msg,const t& state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Restart();

        public:
            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_nat : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_param : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                const Messaging& msg,
                const Reals& reals,
                const Nats& nats,
                const Params& params,
                const X_Vectors& xs,
                const Y_Vectors& ys
            ) {
                peopt::checkLabels <is_real>
                    (msg,reals.first," real name: ");
                peopt::checkLabels <is_nat>
                    (msg,nats.first," natural name: ");
                peopt::checkLabels <is_param>
                    (msg,params.first," paramater name: ");
                peopt::checkLabels <is_x>
                    (msg,xs.first," variable name: ");
                peopt::checkLabels <is_y>
                    (msg,ys.first," equality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string,std::string,std::string>
            {
                std::string operator() (
                    std::string label,
                    std::string val
                ) {

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
                typename State::t& state, 
                X_Vectors& xs,
                Y_Vectors& ys
            ) {
                ys.first.push_back("y");
                ys.second.splice(ys.second.end(),state.y);
                ys.first.push_back("dy");
                ys.second.splice(ys.second.end(),state.dy);
                ys.first.push_back("g_x");
                ys.second.splice(ys.second.end(),state.g_x);
                ys.first.push_back("gpxdxn_p_gx");
                ys.second.splice(ys.second.end(),state.gpxdxn_p_gx);
                ys.first.push_back("gpxdxt");
                ys.second.splice(ys.second.end(),state.gpxdxt);
                xs.first.push_back("dx_n");
                xs.second.splice(xs.second.end(),state.dx_n);
                xs.first.push_back("dx_ncp");
                xs.second.splice(xs.second.end(),state.dx_ncp);
                xs.first.push_back("dx_t");
                xs.second.splice(xs.second.end(),state.dx_t);
                xs.first.push_back("dx_t_uncorrected");
                xs.second.splice(xs.second.end(),state.dx_t_uncorrected);
                xs.first.push_back("dx_tcp_uncorrected");
                xs.second.splice(xs.second.end(),state.dx_tcp_uncorrected);
                xs.first.push_back("H_dxn");
                xs.second.splice(xs.second.end(),state.H_dxn);
                xs.first.push_back("W_gradpHdxn");
                xs.second.splice(xs.second.end(),state.W_gradpHdxn);
                xs.first.push_back("H_dxtuncorrected");
                xs.second.splice(xs.second.end(),state.H_dxtuncorrected);
            }

            // Copy out all the scalar information
            static void stateToScalars(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
            ) { 
                // Copy in all the real numbers
                reals.first.push_back("zeta");
                reals.second.push_back(state.zeta);
                reals.first.push_back("eta0");
                reals.second.push_back(state.eta0);
                reals.first.push_back("rho");
                reals.second.push_back(state.rho);
                reals.first.push_back("rho_old");
                reals.second.push_back(state.rho_old);
                reals.first.push_back("rho_bar");
                reals.second.push_back(state.rho_bar);
                reals.first.push_back("eps_constr");
                reals.second.push_back(state.eps_constr);
                reals.first.push_back("xi_qn");
                reals.second.push_back(state.xi_qn);
                reals.first.push_back("xi_pg");
                reals.second.push_back(state.xi_pg);
                reals.first.push_back("xi_proj");
                reals.second.push_back(state.xi_proj);
                reals.first.push_back("xi_tang");
                reals.second.push_back(state.xi_tang);
                reals.first.push_back("xi_lmh");
                reals.second.push_back(state.xi_lmh);
                reals.first.push_back("xi_lmg");
                reals.second.push_back(state.xi_lmg);
                reals.first.push_back("xi_4");
                reals.second.push_back(state.xi_4);
                reals.first.push_back("rpred");
                reals.second.push_back(state.rpred);
                reals.first.push_back("norm_gxtyp");
                reals.second.push_back(state.norm_gxtyp);
                reals.first.push_back("norm_gpxdxnpgx");
                reals.second.push_back(state.norm_gpxdxnpgx);

                // Copy in all the natural numbers
                nats.first.push_back("augsys_iter_max");
                nats.second.push_back(state.augsys_iter_max);
                nats.first.push_back("augsys_rst_freq");
                nats.second.push_back(state.augsys_rst_freq);

                // Copy in all the parameters
                params.first.push_back("PSchur_left_type");
                params.second.push_back(
                    Operators::to_string(state.PSchur_left_type));
                params.first.push_back("PSchur_right_type");
                params.second.push_back(
                    Operators::to_string(state.PSchur_right_type));
            }
            
            // Copy in all equality multipliers 
            static void vectorsToState(
                typename State::t& state,
                X_Vectors& xs,
                Y_Vectors& ys
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
                        state.y.splice(state.y.end(),ys.second,y0);
                    else if(*name0=="dy")
                        state.dy.splice(state.dy.end(),ys.second,y0);
                    else if(*name0=="g_x")
                        state.g_x.splice(state.g_x.end(),ys.second,y0);
                    else if(*name0=="gpxdxn_p_gx")
                        state.gpxdxn_p_gx.splice(
                            state.gpxdxn_p_gx.end(),ys.second,y0);
                    else if(*name0=="gpxdxt")
                        state.gpxdxt.splice(state.gpxdxt.end(),ys.second,y0);
                    
                    // Remove the string corresponding to the element just
                    // spliced if splicing occured.
                    if(ys.first.size() != ys.second.size())
                        ys.first.erase(name0);
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
                        state.dx_n.splice(state.dx_n.end(),xs.second,x0);
                    else if(*name0=="dx_ncp")
                        state.dx_ncp.splice(state.dx_ncp.end(),xs.second,x0);
                    else if(*name0=="dx_t")
                        state.dx_t.splice(state.dx_t.end(),xs.second,x0);
                    else if(*name0=="dx_t_uncorrected")
                        state.dx_t_uncorrected.splice(
                            state.dx_t_uncorrected.end(),xs.second,x0);
                    else if(*name0=="dx_tcp_uncorrected")
                        state.dx_tcp_uncorrected.splice(
                            state.dx_tcp_uncorrected.end(),xs.second,x0);
                    else if(*name0=="H_dxn")
                        state.H_dxn.splice(
                            state.H_dxn.end(),xs.second,x0);
                    else if(*name0=="W_gradpHdxn")
                        state.W_gradpHdxn.splice(
                            state.W_gradpHdxn.end(),xs.second,x0);
                    else if(*name0=="H_dxtuncorrected")
                        state.H_dxtuncorrected.splice(
                            state.H_dxtuncorrected.end(),xs.second,x0);

                    // Remove the string corresponding to the element just
                    // spliced if splicing occured.
                    if(xs.first.size() != xs.second.size())
                        xs.first.erase(name0);
                }
            }
            
            // Copy in all the scalar information
            static void scalarsToState(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
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
                typename State::t& state,
                X_Vectors& xs,
                Y_Vectors& ys,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::Restart::stateToVectors(state,xs);
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
                const Messaging& msg,
                typename State::t& state,
                X_Vectors& xs,
                Y_Vectors& ys,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,ys);

                // Check the strings used to represent parameters
                checkParams <checkParamVal> (msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::Restart::vectorsToState(state,xs);
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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Functions();

        public:
            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of auto_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                t& operator = (const t&);
                t(const t&);

            public:
                // Equality constraints 
                std::auto_ptr <VectorValuedFunction <Real,XX,YY> > g;

                // Left preconditioner for the augmented system
                std::auto_ptr <Operator <Real,YY,YY> > PSchur_left;

                // Right preconditioner for the augmented system
                std::auto_ptr <Operator <Real,YY,YY> > PSchur_right;
                
                // Initialize all of the pointers to null
                t() : Unconstrained <Real,XX>::Functions::t(), g(NULL),
                    PSchur_left(NULL), PSchur_right(NULL) {}
            };

            struct EqualityModifications
                : public peopt::ScalarValuedFunctionModifications <Real,XX>
            {
            private:
                // Underlying modification.  This takes control of the memory
                std::auto_ptr <
                    peopt::ScalarValuedFunctionModifications <Real,XX> > f_mod;

                // Equality constraint.
                const peopt::VectorValuedFunction <Real,XX,YY>& g;

                // Reference to equality Lagrange multiplier
                const Y_Vector& y;

                // Reference to parameter for the augmented-Lagrangian
                const Real& rho;

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
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_lag
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
                
                // Disallow the default and copy constructors as the assignment
                // operator
                EqualityModifications() {}
                EqualityModifications(const EqualityModifications&);
                EqualityModifications& operator =(const EqualityModifications&);
            public:
                EqualityModifications(
                    const typename State::t& state,
                    typename Functions::t& fns
                ) : f_mod(fns.f_mod),
                    g(*(fns.g)),
                    y(state.y.back()),
                    rho(state.rho)
                { 
                    // Create some shortcuts
                    const X_Vector& x=state.x.back();
                    const Y_Vector& y=state.y.back();

                    // Allocate a bit of memory for the workspace
                    X::init(x,grad_tmp); 
                    X::init(x,x_tmp1); 
                    Y::init(y,y_tmp1); 

                    // Allocate memory for the caching
                    X::init(x,x_merit.second);
                        x_merit.first=false;
                    Y::init(y,g_x);
                    X::init(x,x_grad.second);
                        x_grad.first=false;
                    Y::init(y,y_grad.second);
                        y_grad.first=false;
                    X::init(x,gpxsy);
                }

                // Merit function additions to the objective
                virtual Real merit(const X_Vector& x,const Real& f_x) const {
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
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_stop
                ) const {
                    f_mod->grad_stop(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_stop);
                }

                // Diagnostic modification of the gradient
                virtual void grad_diag(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_diag
                ) const {
                    f_mod->grad_diag(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_diag);
                }

                // Modification of the gradient when finding a trial step
                virtual void grad_step(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_step
                ) const {
                    f_mod->grad_step(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_step);
                }

                // Modification of the gradient for a quasi-Newton method 
                virtual void grad_quasi(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_quasi
                ) const {
                    f_mod->grad_quasi(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_quasi);
                }

                // Modification of the gradient when solving for the equality
                // multiplier
                virtual void grad_mult(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_mult
                ) const {
                    f_mod->grad_mult(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_mult);
                }

                // Modification of the Hessian-vector product when finding a
                // trial step
                virtual void hessvec_step(
                    const X_Vector& x,
                    const X_Vector& dx,
                    const X_Vector& H_dx,
                    X_Vector& Hdx_step 
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
                void operator () (const Y_Vector& dy,Y_Vector& result) const{
                    Y::copy(dy,result);
                }
            };

            // Check that all the functions are defined
            static void check(const Messaging& msg,const t& fns) {

                // Check the unconstrained pieces
                Unconstrained <Real,XX>::Functions::check(msg,fns);
                
                // Check that the equality constraints exist 
                if(fns.g.get()==NULL)
                    msg.error("Missing the equality constraint definition.");

                // Check that preconditioners exist
                if(fns.PSchur_left.get()==NULL)
                    msg.error("Missing a left preconditioner for the "
                        "augmented system.");
                if(fns.PSchur_right.get()==NULL)
                    msg.error("Missing a right preconditioner for the "
                        "augmented system.");
            }

            // Initialize any missing functions for just equality constrained 
            // optimization.
            static void init_(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
            ) {
                // Determine the left preconditioner for the augmented system
                switch(state.PSchur_left_type){
                    case Operators::Identity:
                        fns.PSchur_left.reset(new Identity());
                        break;
                    case Operators::UserDefined:
                        if(fns.PSchur_left.get()==NULL)
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
                        if(fns.PSchur_right.get()==NULL)
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
                if(fns.TRS.get()==NULL)
                    fns.TRS.reset(new typename
                        Unconstrained <Real,XX>::Functions::Identity());
                
                // Check that all functions are defined 
                check(msg,fns);
                
                // Modify the objective 
                fns.f_mod.reset(new EqualityModifications(state,fns));
            }

            // Initialize any missing functions 
            static void init(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
            ) {
                Unconstrained <Real,XX>
                    ::Functions::init_(msg,state,fns);
                EqualityConstrained <Real,XX,YY>
                    ::Functions::init_(msg,state,fns);
            }
        };
        
        // Contains functions that assist in creating an output for diagonstics
        struct Printer {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Printer();

        public:
            // Gets the header for the state information
            static void getStateHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) { 
                // Norm of the constrained 
                out.push_back(atos <> ("||g(x)||"));
                    
                // Trust-region information
                out.push_back(atos <> ("ared"));
                out.push_back(atos <> ("pred"));
                out.push_back(atos <> ("ared/pred"));
                   
                // Krylov method information
                out.push_back(atos <> ("KryIter"));
                out.push_back(atos <> ("KryErr"));
                out.push_back(atos <> ("KryWhy"));
            }
            // Combines all of the state headers
            static void getStateHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer::getStateHeader_
                    (state,out);
            }

            // Gets the state information for output
            static void getState_(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                // Create some shortcuts
                const Y_Vector& g_x = state.g_x.back();
                const Natural& krylov_iter=state.krylov_iter;
                const Real& krylov_rel_err=state.krylov_rel_err;
                const KrylovStop::t& krylov_stop=state.krylov_stop;
                const Natural& iter=state.iter;
                const Natural& rejected_trustregion=state.rejected_trustregion;
                const Real& pred = state.pred;
                const Real& ared = state.ared;

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
                out.push_back(atos <> (norm_gx));
                    
                // Actual vs. predicted reduction 
                if(!opt_begin) {
                    out.push_back(atos <> (ared));
                    out.push_back(atos <> (pred));
                    out.push_back(atos <> (ared/pred));
                } else 
                    for(Natural i=0;i<3;i++) out.push_back(blankSeparator);
                
                // Krylov method information
                if(!opt_begin) {
                    out.push_back(atos <> (krylov_iter));
                    out.push_back(atos <> (krylov_rel_err));
                    out.push_back(atos <> (krylov_stop));
                } else 
                    for(Natural i=0;i<3;i++) out.push_back(blankSeparator);

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getState_(fns,state,blank,out);
            }
            
            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov information
            static void getKrylov(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };

        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Algorithms();

        public:
            // The operator for the augmented system,
            //
            // [ I      g'(x)* ]
            // [ g'(x)  0      ]
            //
            struct AugmentedSystem: public Operator <Real,XXxYY,XXxYY> {
            private:
                const typename State::t& state;
                const typename Functions::t& fns;
                const X_Vector& x_base;
            public:
                AugmentedSystem(
                    const typename State::t& state_,
                    const typename Functions::t& fns_,
                    const X_Vector& x_base_
                ) : state(state_), fns(fns_), x_base(x_base_) {}
                
                // Operator interface
                void operator () (
                    const XxY_Vector& dx_dy,
                    XxY_Vector& result
                ) const{
                    // Create some shortcuts
                    const VectorValuedFunction <Real,XX,YY>& g=*(fns.g);

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
                const Operator <Real,XX,XX>& PH_x;
                const Operator <Real,YY,YY>& PH_y;
            public:
                BlockDiagonalPreconditioner(
                    const Operator <Real,XX,XX>& PH_x_,
                    const Operator <Real,YY,YY>& PH_y_ 
                ) : PH_x(PH_x_), PH_y(PH_y_) {}
                
                // Operator interface
                void operator () (
                    const XxY_Vector& dx_dy,
                    XxY_Vector& result
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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                explicit QNManipulator(
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    const Natural& iter,
                    const typename XXxYY <Real>::Vector& xx,
                    const typename XXxYY <Real>::Vector& bb,
                    Real& eps
                ) const {
                    // Create some shortcuts
                    const Real& xi_qn = state.xi_qn;
                    const Real& norm_gxtyp = state.norm_gxtyp;
                    const Real& eps_constr= state.eps_constr;

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const VectorValuedFunction <Real,XX,YY>& g=*(fns.g);
                const X_Vector& x=state.x.front();
                const Y_Vector& g_x=state.g_x.front();
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                const Real& delta = state.delta;
                const Real& zeta = state.zeta;
                X_Vector& dx_ncp=state.dx_ncp.front();
                X_Vector& dx_n=state.dx_n.front();

                // Find the Cauchy point.

                // Find g'(x)*g(x)
                X_Vector gps_g; X::init(x,gps_g);
                g.ps(x,g_x,gps_g);

                // Find g'(x)g'(x)*g(x)
                Y_Vector gp_gps_g; Y::init(g_x,gp_gps_g);
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
                XxY_Vector x0; X::init(x,x0.first); Y::init(g_x,x0.second);
                XxY::zero(x0);

                // Create the rhs, b0=(-dx_ncp,-g'(x)dx_ncp-g(x)) 
                XxY_Vector b0; XxY::init(x0,b0);
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
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                X_Vector& dx_dnewton = x0.first;

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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                explicit NullspaceProjForGradLagPlusHdxnManipulator(
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    const Natural& iter,
                    const typename XXxYY <Real>::Vector& xx,
                    const typename XXxYY <Real>::Vector& bb,
                    Real& eps
                ) const {
                    // Create some shortcuts
                    const Real& xi_pg = state.xi_pg;
                    const Real& delta = state.delta;

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& dx_n=state.dx_n.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& H_dxn=state.H_dxn.front();
                const Y_Vector& y=state.y.front();
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                X_Vector& W_gradpHdxn=state.W_gradpHdxn.front();

                // Find the gradient modifications for the step computation
                X_Vector grad_step;
                    X::init(grad,grad_step);
                    f_mod.grad_step(x,grad,grad_step);
               
                // Add the Hessian modifications to H(x)dx_n
                X_Vector Hdxn_step;
                    X::init(x,Hdxn_step);
                    f_mod.hessvec_step(x,dx_n,H_dxn,Hdxn_step);

                // grad_p_Hdxn <- H dxn_step
                X_Vector grad_p_Hdxn;
                    X::init(x,grad_p_Hdxn);
                    X::copy(Hdxn_step,grad_p_Hdxn);

                // grad_p_Hdxn <- grad f(x) + H dx_n
                X::axpy(Real(1.),grad_step,grad_p_Hdxn);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0;
                    X::init(x,x0.first);
                    Y::init(y,x0.second);
                    XxY::zero(x0);

                // Create the rhs, b0=(grad f(x) + H dx_n,0)
                XxY_Vector b0;
                    XxY::init(x0,b0);
                    X::copy(grad_p_Hdxn,b0.first);
                    Y::zero(b0.second);
            
                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l (I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r (I,*(fns.PSchur_right));

                // Solve the augmented system for the nullspace projection 
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                explicit NullspaceProjForKrylovMethodManipulator (
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    const Natural& iter,
                    const typename XXxYY <Real>::Vector& xx,
                    const typename XXxYY <Real>::Vector& bb,
                    Real& eps
                ) const {
                    // Create some shortcuts
                    const Real& xi_proj = state.xi_proj;

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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                NullspaceProjForKrylovMethod(
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
               
                // Project dx_t into the nullspace of g'(x)
                void operator () (
                    const X_Vector& dx_t_uncorrected,
                    X_Vector& result
                ) const{
                    // Create some shortcuts
                    const X_Vector& x=state.x.front();
                    const Y_Vector& y=state.y.front();
                    const unsigned int augsys_iter_max=state.augsys_iter_max;
                    const unsigned int augsys_rst_freq=state.augsys_rst_freq;

                    // Create the initial guess, x0=(0,0)
                    XxY_Vector x0;
                        X::init(x,x0.first);
                        Y::init(y,x0.second);
                        XxY::zero(x0);

                    // Create the rhs, b0=(dx_t_uncorrected,0)
                    XxY_Vector b0;
                        XxY::init(x0,b0);
                        X::copy(dx_t_uncorrected,b0.first);
                        Y::zero(b0.second);
                
                    // Build Schur style preconditioners
                    typename Unconstrained <Real,XX>::Functions::Identity I;
                    BlockDiagonalPreconditioner
                        PAugSys_l(I,*(fns.PSchur_left));
                    BlockDiagonalPreconditioner
                        PAugSys_r(I,*(fns.PSchur_right));

                    // Solve the augmented system for the nullspace projection 
                    std::pair <Real,Natural> err_iter
                        = peopt::gmres <Real,XXxYY> (
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& dx_n=state.dx_n.front();
                const X_Vector& W_gradpHdxn=state.W_gradpHdxn.front();
                const Real& delta = state.delta;
                const Real& eps_krylov=state.eps_krylov;
                const Natural& krylov_iter_max=state.krylov_iter_max;
                const Natural& krylov_orthog_max=state.krylov_orthog_max;
                const KrylovSolverTruncated::t& krylov_solver
                    = state.krylov_solver;
                X_Vector& dx_t_uncorrected=state.dx_t_uncorrected.front();
                X_Vector& dx_tcp_uncorrected=state.dx_tcp_uncorrected.front();
                Real& krylov_rel_err=state.krylov_rel_err;
                Natural& krylov_iter=state.krylov_iter;
                Natural& krylov_iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const Operator <Real,XX,XX>& TRS=*(fns.TRS);
                    
                // Setup the Hessian operator and allocate memory for the
                // Cauchy point.
                typename Unconstrained <Real,XX>::Algorithms::HessianOperator
                    H(f,f_mod,x);

                // Find the quantity - W (g + H dxn).  We use this as the
                // RHS in the linear system solve.
                X_Vector minus_W_gradpHdxn;
                    X::init(x,minus_W_gradpHdxn);
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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                TangentialStepManipulator (
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    const Natural& iter,
                    const typename XXxYY <Real>::Vector& xx,
                    const typename XXxYY <Real>::Vector& bb,
                    Real& eps
                ) const {
                    // Create some shortcuts
                    const X_Vector& dx_n=state.dx_n.front();
                    const Real& xi_tang = state.xi_tang;
                    const Real& delta = state.delta;

                    // dxn_p_dxt <- dx_n + dx_t
                    X_Vector dxn_p_dxt; X::init(dx_n,dxn_p_dxt);
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const X_Vector& x=state.x.front();
                const Y_Vector& y=state.y.front();
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                const X_Vector& dx_t_uncorrected=state.dx_t_uncorrected.front();
                X_Vector& dx_t=state.dx_t.front();

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0;
                    X::init(x,x0.first);
                    Y::init(y,x0.second);
                    XxY::zero(x0);

                // Create the rhs, b0=(dx_t_uncorrected,0);
                XxY_Vector b0;
                    XxY::init(x0,b0);
                    X::copy(dx_t_uncorrected,b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the tangential step 
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                const typename State::t& state;
                const typename Functions::t& fns;
            public:
                LagrangeMultiplierStepManipulator (
                    const typename State::t& state_,
                    const typename Functions::t& fns_
                ) : state(state_), fns(fns_) {}
                void operator () (
                    const Natural& iter,
                    const typename XXxYY <Real>::Vector& xx,
                    const typename XXxYY <Real>::Vector& bb,
                    Real& eps
                ) const {
                    // Create some shortcuts
                    const Real& xi_lmh = state.xi_lmh;
                    const Real& xi_lmg = state.xi_lmg;
                
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                const X_Vector& grad=state.grad.front();
                Y_Vector& y=state.y.front();

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_mult;
                    X::init(grad,grad_mult);
                    f_mod.grad_mult(x,grad,grad_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0;
                    X::init(x,x0.first);
                    Y::init(y,x0.second);
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x,y),0);
                XxY_Vector b0;
                    XxY::init(x0,b0);
                    X::copy(grad_mult,b0.first);
                    X::scal(Real(-1.),b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the initial Lagrange
                // multiplier 
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& grad=state.grad.front();
                const X_Vector& dx=state.dx.front();
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                X_Vector& x=state.x.front();
                Y_Vector& dy=state.dy.front();

                // x_p_dx <- x + dx
                X_Vector x_p_dx;
                    X::init(x,x_p_dx);
                    X::copy(x,x_p_dx);
                    X::axpy(Real(1.),dx,x_p_dx);

                // grad_xpdx <- L(x+dx,y) = grad f(x+dx) + g'(x+dx)*y
                X_Vector grad_xpdx; 
                    X::init(x,grad_xpdx);
                    f.grad(x_p_dx,grad_xpdx);

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_xpdx_mult;
                    X::init(grad,grad_xpdx_mult);
                    f_mod.grad_mult(x_p_dx,grad_xpdx,grad_xpdx_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0;
                    X::init(x,x0.first);
                    Y::init(dy,x0.second);
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x+dx,y),0);
                XxY_Vector b0;
                    XxY::init(x0,b0);
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
                X_Vector x_save;
                    X::init(x,x_save);
                    X::copy(x,x_save);
                X::copy(x_p_dx,x);

                // Solve the augmented system for the Lagrange multiplier step 
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& grad=state.grad.front(); 
                const Natural& augsys_iter_max=state.augsys_iter_max;
                const Natural& augsys_rst_freq=state.augsys_rst_freq;
                X_Vector& x=state.x.front();
                Y_Vector& dy=state.dy.front();

                // grad_x <- L(x,y) = grad f(x) + g'(x)*y
                X_Vector grad_x; 
                    X::init(x,grad_x);
                    f.grad(x,grad_x);

                // Find the gradient modifications for the Lagrange multiplier
                // computation
                X_Vector grad_x_mult;
                    X::init(grad,grad_x_mult);
                    f_mod.grad_mult(x,grad_x,grad_x_mult);

                // Create the initial guess, x0=(0,0)
                XxY_Vector x0;
                    X::init(x,x0.first);
                    Y::init(dy,x0.second);
                    XxY::zero(x0);

                // Create the rhs, b0=(-grad L(x,y),0);
                XxY_Vector b0;
                    XxY::init(x0,b0);
                    X::copy(grad_x_mult,b0.first);
                    X::scal(Real(-1.),b0.first);
                    Y::zero(b0.second);

                // Build Schur style preconditioners
                typename Unconstrained <Real,XX>::Functions::Identity I;
                BlockDiagonalPreconditioner PAugSys_l(I,*(fns.PSchur_left));
                BlockDiagonalPreconditioner PAugSys_r(I,*(fns.PSchur_right));

                // Solve the augmented system for the Lagrange multiplier step 
                std::pair <Real,Natural> err_iter = peopt::gmres <Real,XXxYY> (
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& grad=state.grad.front();
                const X_Vector& W_gradpHdxn=state.W_gradpHdxn.front();
                const X_Vector& H_dxn=state.H_dxn.front();
                const X_Vector& dx_n=state.dx_n.front();
                const X_Vector& H_dxtuncorrected=state.H_dxtuncorrected.front();
                const X_Vector& dx_t_uncorrected=state.dx_t_uncorrected.front();
                const Y_Vector& gpxdxn_p_gx=state.gpxdxn_p_gx.front();
                const Y_Vector& dy=state.dy.front();
                const Y_Vector& g_x = state.g_x.front();
                const Real& rho=state.rho;
                const Real& norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                Real& pred=state.pred;
                
                // Find || g(x) ||
                Real norm_gx = sqrt(Y::innr(g_x,g_x));
                
                // Find the gradient modifications for step computation 
                X_Vector grad_step;
                    X::init(grad,grad_step);
                    f_mod.grad_step(x,grad,grad_step);
                
                // Add the Hessian modifications to H(x)dx_n
                X_Vector Hdxn_step;
                    X::init(x,Hdxn_step);
                    f_mod.hessvec_step(x,dx_n,H_dxn,Hdxn_step);
                
                // Add the Hessian modifications to H(x)dx_t_uncorrected
                X_Vector H_dxtuncorrected_step;
                    X::init(x,H_dxtuncorrected_step);
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Y_Vector& g_x=state.g_x.back();
                const Real& pred=state.pred;
                const Real& norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                const Real& rho_old=state.rho_old;
                const Real& rho_bar=state.rho_bar;
                Real& rho=state.rho;
               
                // norm_gx <- || g(x) ||
                const Real& norm_gx=sqrt(Y::innr(g_x,g_x));

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Y_Vector& dy=state.dy.front();
                const Y_Vector& gpxdxn_p_gx=state.gpxdxn_p_gx.front();
                const Y_Vector& gpxdxt=state.gpxdxt.front();
                const Real& rho=state.rho;
                Real& rpred=state.rpred;

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
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create shortcuts to some elements in the state
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.front();
                const X_Vector& dx=state.dx.front();
                const Y_Vector& dy=state.dy.front();
                const Real& eta1=state.eta1;
                const Real& eta2=state.eta2;
                const Real& f_x=state.f_x;
                const KrylovStop::t& krylov_stop=state.krylov_stop;
                Y_Vector& y=state.y.front();
                Real& delta=state.delta;
                Real& ared=state.ared;
                Real& pred=state.pred;
                Real& f_xpdx=state.f_xpdx;
                
                // Allocate memory for temporaries that we need
                X_Vector x_p_dx; X::init(x,x_p_dx);

                // Determine x+dx 
                X::copy(dx,x_p_dx);
                X::axpy(Real(1.),x,x_p_dx);

                // Save the old Lagrange multiplier
                Y_Vector y_old;
                    Y::init(y,y_old);
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
                const Messaging& msg,
                const StateManipulator <EqualityConstrained <Real,XX,YY> >&
                    smanip,
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create some shortcuts
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const VectorValuedFunction <Real,XX,YY>& g=*(fns.g);
                const X_Vector& x=state.x.front();
                const X_Vector& dx_n=state.dx_n.front();
                const X_Vector& dx_t=state.dx_t.front();
                const X_Vector& dx_tcp_uncorrected
                    =state.dx_tcp_uncorrected.front();
                const Real& xi_4=state.xi_4;
                const Real& eta0=state.eta0;
                const Real& eps_dx=state.eps_dx;
                const Real& norm_dxtyp=state.norm_dxtyp;
                const Real& rho_old=state.rho_old;
                const Natural& history_reset=state.history_reset;
                X_Vector& dx=state.dx.front();
                X_Vector& dx_t_uncorrected=state.dx_t_uncorrected.front();
                X_Vector& H_dxn=state.H_dxn.front();
                X_Vector& H_dxtuncorrected=state.H_dxtuncorrected.front();
                Y_Vector& g_x=state.g_x.front();
                Y_Vector& gpxdxn_p_gx=state.gpxdxn_p_gx.front();
                Y_Vector& gpxdxt=state.gpxdxt.front();
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                Real& norm_gpxdxnpgx=state.norm_gpxdxnpgx;
                Real& xi_qn=state.xi_qn;
                Real& xi_pg=state.xi_pg;
                Real& xi_proj=state.xi_proj;
                Real& xi_tang=state.xi_tang;
                Real& xi_lmh=state.xi_lmh;
                Real& pred=state.pred;
                Real& rpred=state.rpred;
                Real& rho=state.rho;
                Real& alpha=state.alpha;
                Natural& rejected_trustregion=state.rejected_trustregion;

                // Create a single temporary vector
                X_Vector x_tmp1; X::init(x,x_tmp1);

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
                    // linesearch parameter, alpha.
                    alpha = Real(1.);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::BeforeActualVersusPredicted);

                    // If need be, shorten the step
                    X::scal(alpha,dx);

                    // If we shorten our step, update our Lagrange multiplier
                    // step
                    if(alpha < Real(1.))
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Y_Vector& g_x = state.g_x.back();
                const Real& eps_constr=state.eps_constr;
                const Real& norm_gxtyp=state.norm_gxtyp; 
                StoppingCondition::t& opt_stop=state.opt_stop;
                
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
                const StateManipulator<ProblemClass>& smanip;

                // A reference to the messaging object
                const Messaging& msg;

            public:
                CompositeStepManipulator(
                    const StateManipulator <ProblemClass>& smanip_,
                    const Messaging& msg_
                ) : smanip(smanip_), msg(msg_) {}


                // Application
                void operator () (
                    const typename ProblemClass::Functions::t& fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t loc
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
                    const typename Functions::t& fns
                        =dynamic_cast <const typename Functions::t&> (fns_);
                    typename State::t& state 
                        =dynamic_cast <typename State::t&> (state_);
                
                    // Create some shortcuts
                    const ScalarValuedFunctionModifications <Real,XX>& f_mod
                        = *(fns.f_mod);
                    const VectorValuedFunction <Real,XX,YY>& g=*(fns.g);
                    const X_Vector& x=state.x.front();
                    const Y_Vector& dy=state.dy.front();
                    const Real& rho = state.rho;
                    const Natural& krylov_iter_max = state.krylov_iter_max;
                    AlgorithmClass::t& algorithm_class=state.algorithm_class;
                    X_Vector& grad=state.grad.front();
                    Y_Vector& y=state.y.front();
                    Y_Vector& g_x=state.g_x.front();
                    Real& norm_gxtyp = state.norm_gxtyp;
                    Real& rho_old = state.rho_old;
                    Real& norm_gradtyp = state.norm_gradtyp;
                    Natural& krylov_orthog_max = state.krylov_orthog_max;

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
                        X_Vector grad_stop;
                            X::init(grad,grad_stop);
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
                const Messaging& msg,
                typename Functions::t& fns,
                typename State::t& state
            ){
                // Create an empty state manipulator
                StateManipulator <EqualityConstrained <Real,XX,YY> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                const Messaging& msg,
                const StateManipulator <EqualityConstrained <Real,XX,YY> >&
                    smanip,
                typename Functions::t& fns,
                typename State::t& state
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
    struct InequalityConstrained : public virtual Unconstrained <Real,XX> {
    private:
        // This is a templated namespace.  Do not allow construction.
        InequalityConstrained();

    public:
        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;
        
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

        // Functions that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            State();

        public:
            // The actual internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                t& operator = (const t&);
                t(const t&);

            public:
                // Lagrange multiplier (dual variable) for the inequality
                // constraints 
                std::list <Z_Vector> z;
                
                // Step in the Lagrange multiplier 
                std::list <Z_Vector> dz;

                // The inequality constraint evaluated at x.  In theory,
                // we can always just evaluate this when we need it.  However,
                // we do require its computation both in the gradient as well
                // as Hessian calculations.  More specifically, when computing
                // with SDP constraints, we require the Schur factorization
                // of this quantity in order to solve the Sylvester equations
                // for Linv(h(x)).  By caching it, we have the ability
                // to cache the Schur factorization.
                std::list <Z_Vector> h_x;

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
                t() {
                    InequalityConstrained <Real,XX,ZZ>::State
                        ::init_params(*this);
                }
                t(const X_Vector& x,const Z_Vector& z) {
                    InequalityConstrained <Real,XX,ZZ>::State
                        ::init_params(*this);
                    InequalityConstrained <Real,XX,ZZ>::State
                        ::init_vectors(*this,x,z);
                }
            };
                
            // This initializes all the parameters required for inequality
            // constrained optimization.  
            static void init_params_(t& state) {
                state.mu = std::numeric_limits<Real>::quiet_NaN();
                state.mu_est = std::numeric_limits<Real>::quiet_NaN();
                state.mu_typ = std::numeric_limits<Real>::quiet_NaN();
                state.eps_mu= Real(1e-8);
                state.sigma = Real(0.5);
                state.gamma = Real(0.95);
                state.ipm = InteriorPointMethod::PrimalDual;
                state.cstrat = CentralityStrategy::Constant; 
            }
            static void init_params(t& state) {
                Unconstrained <Real,XX>::State::init_params_(state); 
                InequalityConstrained <Real,XX,ZZ>::State::init_params_(state);
            }

            // This initializes all the variables required for inequality
            // constrained optimization.  
            static void init_vectors_(
                t& state,
                const X_Vector& x,
                const Z_Vector& z
            ) {
                // Allocate memory for z
                state.z.clear();
                    state.z.push_back(Z_Vector());
                    Z::init(z,state.z.back());
                    Z::copy(z,state.z.back());

                // Allocate memory for dz
                state.dz.clear();
                    state.dz.push_back(Z_Vector());
                    Z::init(z,state.dz.back());

                // Allocate memory for h(x)
                state.h_x.clear();
                    state.h_x.push_back(Z_Vector());
                    Z::init(z,state.h_x.back());
            }
            static void init_vectors(
                t& state,
                const X_Vector& x,
                const Z_Vector& z
            ) {
                Unconstrained <Real,XX>::State::init_vectors_(state,x); 
                InequalityConstrained <Real,XX,ZZ>::State
                    ::init_vectors_(state,x,z); 
            }
           
            // Initializes everything
            static void init(t& state, const X_Vector& x, const Z_Vector& z) {
                init_params(state);
                init_vectors(state,x,z);
            }

            // Check that we have a valid set of parameters.  
            static void check_(const Messaging& msg,const t& state) {
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
            static void check(const Messaging& msg,const t& state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                InequalityConstrained <Real,XX,ZZ>::State::check_(msg,state);
            }
        };
        // Utilities for restarting the optimization
        struct Restart {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Restart();

        public:
            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_nat : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                        ::is_nat()(name)
                    )
                        return true;
                    else
                        return false;
                }
            };
           
            // Checks whether we have a valid parameter label.
            struct is_param : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                const Messaging& msg,
                const Reals& reals,
                const Nats& nats,
                const Params& params,
                const X_Vectors& xs,
                const Z_Vectors& zs
            ) {
                peopt::checkLabels <is_real>
                    (msg,reals.first," real name: ");
                peopt::checkLabels <is_nat>
                    (msg,nats.first," natural name: ");
                peopt::checkLabels <is_param>
                    (msg,params.first," paramater name: ");
                peopt::checkLabels <is_x>
                    (msg,xs.first," variable name: ");
                peopt::checkLabels <is_z>
                    (msg,zs.first," inequality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string,std::string,std::string>
            {
                std::string operator() (
                    std::string label,
                    std::string val
                ) {

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
                typename State::t& state, 
                X_Vectors& xs,
                Z_Vectors& zs
            ) {
                zs.first.push_back("z");
                zs.second.splice(zs.second.end(),state.z);
                zs.first.push_back("dz");
                zs.second.splice(zs.second.end(),state.dz);
                zs.first.push_back("h_x");
                zs.second.splice(zs.second.end(),state.h_x);
            }
            
            // Copy out the scalar information
            static void stateToScalars(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy in all the real numbers
                reals.first.push_back("mu");
                reals.second.push_back(state.mu);
                reals.first.push_back("mu_est");
                reals.second.push_back(state.mu_est);
                reals.first.push_back("mu_typ");
                reals.second.push_back(state.mu_typ);
                reals.first.push_back("eps_mu");
                reals.second.push_back(state.eps_mu);
                reals.first.push_back("sigma");
                reals.second.push_back(state.sigma);
                reals.first.push_back("gamma");
                reals.second.push_back(state.gamma);

                // Copy in all of the parameters
                params.first.push_back("ipm");
                params.second.push_back(
                    InteriorPointMethod::to_string(state.ipm));
                params.first.push_back("cstrat");
                params.second.push_back(
                    CentralityStrategy::to_string(state.cstrat));
            }
            
            // Copy in inequality multipliers 
            static void vectorsToState(
                typename State::t& state,
                X_Vectors& xs,
                Z_Vectors& zs
            ) {
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

                    // Remove the string corresponding to the element just
                    // spliced if splicing occured.
                    if(xs.first.size() != xs.second.size())
                        xs.first.erase(name0);
                }

                typename std::list <X_Vector>::iterator z
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
                        state.z.splice(state.z.end(),zs.second,z0);
                    else if(*name0=="h_x")
                        state.h_x.splice(state.h_x.end(),zs.second,z0);
                    else if(*name0=="dz")
                        state.dz.splice(state.dz.end(),zs.second,z0);

                    // Remove the string corresponding to the element just
                    // spliced if slicing occured.
                    if(zs.first.size() != zs.second.size())
                        zs.first.erase(name0);
                }
            }
            
            // Copy in the scalar information
            static void scalarsToState(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
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
                typename State::t& state,
                X_Vectors& xs,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::Restart::stateToVectors(state,xs);
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
                const Messaging& msg,
                typename State::t& state,
                X_Vectors& xs,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,zs);

                // Check the strings used to represent parameters
                checkParams <checkParamVal> (msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::Restart::vectorsToState(state,xs);
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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Functions();

        public:
            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of auto_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                t& operator = (const t&);
                t(const t&);

            public:
                // Inequality constraints 
                std::auto_ptr <VectorValuedFunction <Real,XX,ZZ> > h;
                
                // Initialize all of the pointers to null
                t() : Unconstrained <Real,XX>::Functions::t(), h(NULL) {}
            };

            struct InequalityModifications
                : public peopt::ScalarValuedFunctionModifications <Real,XX>
            {
            private:
                // Underlying modification.  This takes control of the memory
                std::auto_ptr <
                    peopt::ScalarValuedFunctionModifications <Real,XX> > f_mod;

                // Inequality constraint.
                const peopt::VectorValuedFunction <Real,XX,ZZ>& h;
                
                // Inequality Lagrange multiplier
                const Z_Vector& z;

                // Interior point parameter
                const Real& mu;

                // Inequality constraint evaluated at x
                const Z_Vector& h_x;
                
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
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_lag
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
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_schur
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
                
                // Disallow the default and copy constructors as the assignment
                // operator
                InequalityModifications() {}
                InequalityModifications(const InequalityModifications&);
                InequalityModifications&
                    operator = (const InequalityModifications&);
            public:
                InequalityModifications(
                    const typename State::t& state,
                    typename Functions::t& fns
                ) : f_mod(fns.f_mod),
                    h(*(fns.h)),
                    z(state.z.front()),
                    mu(state.mu),
                    h_x(state.h_x.front())
                { 
                    // Create some shortcuts
                    const X_Vector& x=state.x.back();

                    // Allocate a bit of memory for the workspace 
                    X::init(x,grad_tmp);
                    X::init(x,hess_mod);
                    X::init(x,x_tmp1);
                    Z::init(z,z_tmp1);
                    Z::init(z,z_tmp2);
                    
                    // Allocate memory for the caching
                    X::init(x,x_merit.second);
                        x_merit.first=false;
                    Z::init(z,hx_merit);
                    X::init(x,x_lag.second);
                        x_lag.first=false;
                    Z::init(z,z_lag.second);
                        z_lag.first=false;
                    X::init(x,x_schur.second);
                        x_schur.first=false;
                    Z::init(z,z_schur.second);
                        z_schur.first=false;
                    X::init(x,hpxsz);
                    X::init(x,hpxs_invLhx_e);
                }

                // Merit function additions to the objective
                virtual Real merit(const X_Vector& x,const Real& f_x) const {
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
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_stop
                ) const {
                    f_mod->grad_stop(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_stop);
                }

                // Diagnostic modification of the gradient
                virtual void grad_diag(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_diag
                ) const {
                    f_mod->grad_diag(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_diag);
                }

                // Modification of the gradient when finding a trial step
                virtual void grad_step(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_step
                ) const {
                    f_mod->grad_step(x,grad,grad_tmp);
                    grad_schur(x,grad_tmp,grad_step);
                }

                // Modification of the gradient for a quasi-Newton method 
                virtual void grad_quasi(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_quasi
                ) const {
                    f_mod->grad_quasi(x,grad,grad_quasi);
                }

                // Modification of the gradient when solving for the equality
                // multiplier
                virtual void grad_mult(
                    const X_Vector& x,
                    const X_Vector& grad,
                    X_Vector& grad_mult
                ) const {
                    f_mod->grad_mult(x,grad,grad_tmp);
                    grad_lag(x,grad_tmp,grad_mult);
                }

                // Modification of the Hessian-vector product when finding a
                // trial step
                virtual void hessvec_step(
                    const X_Vector& x,
                    const X_Vector& dx,
                    const X_Vector& H_dx,
                    X_Vector& Hdx_step 
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
            static void check(const Messaging& msg,const t& fns) {

                // Check the unconstrained pieces
                Unconstrained <Real,XX>::Functions::check(msg,fns);
                
                // Check that the inequality constraints exist 
                if(fns.h.get()==NULL)
                    msg.error("Missing the inequality constraint definition.");
            }

            // Initialize any missing functions for just inequality constrained 
            // optimization.
            static void init_(
                const Messaging& msg,
                typename State::t& state,
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
                if(fns.TRS.get()==NULL)
                    fns.TRS.reset(new typename Unconstrained <Real,XX>
                        ::Functions::Identity());
#endif
            }

            // Initialize any missing functions 
            static void init(
                const Messaging& msg,
                typename State::t& state,
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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Printer();

        public:
            // Gets the header for the state information
            static void getStateHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                // Print out the current interior point parameter and
                // the estimate of the interior point parameter.
                out.push_back(atos <> ("mu"));
                out.push_back(atos <> ("mu_est"));
            }

            // Combines all of the state headers
            static void getStateHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer::getStateHeader_
                    (state,out);
            }

            // Gets the state information for output
            static void getState_(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {

                // Create some shortcuts
                const Real& mu=state.mu; 
                const Real& mu_est=state.mu_est; 

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Interior point information
                out.push_back(atos <> (mu));
                out.push_back(atos <> (mu_est));

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)=blankSeparator;
            }

            // Combines all of the state information
            static void getState(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getState_(fns,state,blank,out);
            }
            
            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov information
            static void getKrylov(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };

        // This contains the different algorithms used for optimization 
        struct Algorithms {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Algorithms();

        public:
            // An operator to reshape the trust-region radius and, hopefully,
            // keep us away from the boundary.
            struct TrustRegionScaling : public Operator <Real,XX,XX> {
            private:
                // The function h 
                const VectorValuedFunction <Real,XX,ZZ>& h;

                // The value h(x) 
                const Z_Vector& h_x;

                // The current iterate
                const X_Vector& x;

                // Work vectors
                mutable Z_Vector z_tmp1;
                mutable Z_Vector z_tmp2;

            public:
                TrustRegionScaling(
                    const typename Functions::t& fns,
                    const typename State::t& state
                ) : h(*(fns.h)), h_x(state.h_x.back()), x(state.x.back()) {
                    Z::init(h_x,z_tmp1); 
                    Z::init(h_x,z_tmp2);    
                }

                void operator () (const X_Vector& dx,X_Vector& result) const{
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& h_x=state.h_x.front();
                const X_Vector& x=state.x.front();
                const X_Vector& dx=state.dx.front();
                const Real& mu=state.mu;
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                Z_Vector& z=state.z.front();

                // z_tmp1 <- h'(x)dx
                Z_Vector z_tmp1;
                    Z::init(z,z_tmp1);
                    h.p(x,dx,z_tmp1);

                // z_tmp2 <- h'(x)dx o z
                Z_Vector z_tmp2;
                    Z::init(z,z_tmp2);
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& h_x=state.h_x.front();
                const Real& mu=state.mu;
                Z_Vector& z=state.z.front();

                // z_tmp1 <- e 
                Z_Vector z_tmp1;
                    Z::init(z,z_tmp1);
                    Z::id(z_tmp1);

                // z <- inv(L(h(x))) e 
                Z::linv(h_x,z_tmp1,z);

                // z <- mu inv(L(h(x))) e 
                Z::scal(mu,z);
                
                // Symmetrize the iterate 
                Z::symm(z);
            }

#if 1
            // Finds the new inequality Lagrange multiplier step
            // dz = -z + inv L(h(x)) (-h'(x)dx o z + mu e)
            static void findInequalityMultiplierStep(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();
                const X_Vector& x=state.x.front();
                const X_Vector& dx=state.dx.front();
                const Real& mu=state.mu;
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                Z_Vector& dz=state.dz.front();

                // z_tmp1 <- h'(x)dx
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h.p(x,dx,z_tmp1);

                // z_tmp2 <- h'(x)dx o z
                Z_Vector z_tmp2; Z::init(z,z_tmp2);
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
#else
            // Finds the new inequality Lagrange multiplier step
            // dz = -z + inv L(h(x)) (-h'(x)dx o z + mu e)
            static void findInequalityMultiplierStep(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();
                const X_Vector& x=state.x.front();
                const X_Vector& dx=state.dx.front();
                const Real& mu=state.mu;
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                Z_Vector& dz=state.dz.front();

                // z_tmp1 <- h'(x)dx
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h.p(x,dx,z_tmp1);
                
                // z_tmp1 <- -0.5 h'(x)dx
                Z::scal(Real(-0.5),z_tmp1);

                // z_tmp2 <- -0.5 h'(x)dx o z
                Z_Vector z_tmp2; Z::init(z,z_tmp2);
                Z::prod(z_tmp1,z,z_tmp2);

                // dz <- e
                Z::id(dz);

                // z_tmp2 <- -0.5 h'(x)dx o z + mu e
                Z::axpy(mu,dz,z_tmp2);

                // dz <- inv L(h(x)) (-0.5 h'(x)dx o z + mu e)
                Z::linv(h_x,z_tmp2,dz);
               
                // z_tmp3 <- e
                Z_Vector z_tmp3; Z::init(z,z_tmp3);
                Z::id(z_tmp3);

                // z_tmp2 <- inv(L(h(x))) e
                Z::linv(h_x,z_tmp3,z_tmp2);

                // z_tmp3 <- -0.5 h'(x)dx o inv(L(h(x))) e
                Z::prod(z_tmp1,z_tmp2,z_tmp3);

                // z_tmp1 <- z o (-0.5 h'(x)dx o inv(L(h(x))) e)
                Z::prod(z,z_tmp3,z_tmp1);
                
                // dz <- inv L(h(x)) (-0.5 h'(x)dx o z + mu e)
                //     + z o (-0.5 h'(x)dx o inv(L(h(x))) e)
                Z::axpy(Real(1.),z_tmp1,dz);

                // dz <- -z + inv L(h(x)) (-h'(x)dx o z + mu e)
                //          + z o (-0.5 h'(x)dx o inv(L(h(x))) e)
                Z::axpy(Real(-1.),z,dz);
            }
#endif

            // Estimates the interior point parameter with the formula
            // mu = <z,h(x)>/m
            static void estimateInteriorPointParameter(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();
                Real& mu_est=state.mu_est;

                // Determine the scaling factor for the interior-
                // point parameter estimate
                Z_Vector z_tmp; Z::init(z,z_tmp);
                Z::id(z_tmp);
                Real m = Z::innr(z_tmp,z_tmp);

                // Estimate the interior-point parameter
                mu_est = Z::innr(z,h_x) / m;
            }

            // Find interior point parameter
            static void findInteriorPointParameter(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const ScalarValuedFunctionModifications <Real,XX>& f_mod
                    = *(fns.f_mod);
                const X_Vector& x=state.x.back();
                const X_Vector& grad=state.grad.back();
                const Real& mu_est=state.mu_est;
                const Real& mu_typ=state.mu_typ;
                const Real& sigma=state.sigma;
                const Real& eps_mu=state.eps_mu;
                const CentralityStrategy::t& cstrat=state.cstrat;
                const Real& norm_gradtyp=state.norm_gradtyp;
                const Real& eps_grad=state.eps_grad;
                const Natural& iter=state.iter;
                const Real& f_x=state.f_x;
                Real& mu=state.mu;
               
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
                    X_Vector grad_stop;
                        X::init(grad,grad_stop);
                        f_mod.grad_stop(x,grad,grad_stop);
                    const Real norm_grad=sqrt(X::innr(grad_stop,grad_stop));

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Real& mu_est=state.mu_est;
                const Real& mu_typ=state.mu_typ;
                const Real& eps_mu=state.eps_mu;
                const CentralityStrategy::t& cstrat=state.cstrat;
                const Natural& iter=state.iter;
                StoppingCondition::t& opt_stop=state.opt_stop;

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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const Real& gamma=state.gamma;
                const AlgorithmClass::t& algorithm_class
                    =state.algorithm_class;
                const Z_Vector& z=state.z.front();
                const X_Vector& x=state.x.front();
                const Z_Vector& h_x=state.h_x.front();
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                X_Vector& dx=state.dx.front();
                Z_Vector& dz=state.dz.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_; X::init(x,dx_);
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1;
                    X::init(x,x_tmp1);
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1;
                    Z::init(z,z_tmp1);
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
                    alpha *= beta_x;
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
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const Real& gamma=state.gamma;
                const Real& mu=state.mu;
                const AlgorithmClass::t& algorithm_class =state.algorithm_class;
                const Z_Vector& z=state.z.front();
                const X_Vector& x=state.x.front();
                const Z_Vector& h_x=state.h_x.front();
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                X_Vector& dx=state.dx.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_; X::init(x,dx_);
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1;
                    X::init(x,x_tmp1);
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
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
                Z_Vector z_tmp2;
                    Z::init(z,z_tmp2);
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
                    alpha *= alpha0;

            }
            // Conduct a line search that preserves positivity of the
            // primal variable. 
            static void positivityLineSearchLogBarrier(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts 
                const Real& gamma=state.gamma;
                const AlgorithmClass::t& algorithm_class
                    =state.algorithm_class;
                const X_Vector& x=state.x.front();
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                X_Vector& dx=state.dx.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector dx_; X::init(x,dx_);
                X::copy(dx,dx_);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,dx_);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+dx
                X_Vector x_tmp1;
                    X::init(x,x_tmp1);
                    X::copy(x,x_tmp1);
                    X::axpy(Real(1.),dx_,x_tmp1);

                // z_tmp1=h(x+dx)
                Z_Vector z_tmp1;
                    Z::init(z,z_tmp1);
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
                    alpha *= beta_x;
            }

            // This adds the interior point through use of a state manipulator.
            template <typename ProblemClass>
            struct InteriorPointManipulator
                : public StateManipulator <ProblemClass>
            {
            private:
                // A reference to the user-defined state manipulator
                const StateManipulator<ProblemClass>& smanip;

            public:
                InteriorPointManipulator(
                    const StateManipulator <ProblemClass>& smanip_
                ) : smanip(smanip_) {}


                // Application
                void operator () (
                    const typename ProblemClass::Functions::t& fns_,
                    typename ProblemClass::State::t& state_,
                    OptimizationLocation::t loc
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
                    const typename Functions::t& fns
                        =dynamic_cast <const typename Functions::t&> (fns_);
                    typename State::t& state 
                        =dynamic_cast <typename State::t&> (state_);

                    // Create some shorcuts
                    const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                    const X_Vector& x=state.x.front();
                    const InteriorPointMethod::t& ipm=state.ipm;
                    const Real& mu_est = state.mu_est;
                    Z_Vector& z=state.z.front();
                    Z_Vector& h_x=state.h_x.front();
                    Z_Vector& dz=state.dz.front();
                    Real& mu_typ = state.mu_typ;

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
                const Messaging& msg,
                typename Functions::t& fns,
                typename State::t& state
            ){
                // Create an empty state manipulator
                StateManipulator <InequalityConstrained <Real,XX,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                const Messaging& msg,
                const StateManipulator <InequalityConstrained <Real,XX,ZZ> >&
                    smanip,
                typename Functions::t& fns,
                typename State::t& state
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
    private:
        // This is a templated namespace.  Do not allow construction.
        Constrained();

    public:
        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;
        
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

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            State();

        public:
            // The actual internal state of the optimization
            struct t: 
                public EqualityConstrained <Real,XX,YY>::State::t,
                public InequalityConstrained <Real,XX,ZZ>::State::t
            {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Basically, the state can hold a large amount
                // of memory and the safe way to move this memory around
                // is through the use of the capture and release methodology
                // inside of the restart section.
                t& operator = (const t&);
                t(const t&);

            public:
                // Initialization constructors
                t() {
                    Constrained <Real,XX,YY,ZZ>::State::init_params(*this);
                }
                explicit t(
                    const X_Vector& x,const Y_Vector& y,const Z_Vector& z
                ) {
                    Constrained <Real,XX,YY,ZZ>::State::init_params(*this);
                    Constrained <Real,XX,YY,ZZ>::State
                        ::init_vectors(*this,x,y,z);
                }
            };
            
            // This initializes all the parameters required for constrained
            // optimization.  
            static void init_params(t& state) {
                Unconstrained <Real,XX>::State::init_params_(state); 
                EqualityConstrained <Real,XX,YY>::State::init_params_(state);
                InequalityConstrained <Real,XX,ZZ>::State::init_params_(state);
            }

            // This initializes all the variables required for inequality
            // constrained optimization.  
            static void init_vectors(
                t& state,
                const X_Vector& x,
                const Y_Vector& y,
                const Z_Vector& z
            ) {
                Unconstrained <Real,XX>::State::init_vectors_(state,x); 
                EqualityConstrained <Real,XX,YY>::State
                    ::init_vectors_(state,x,y);
                InequalityConstrained <Real,XX,ZZ>::State
                    ::init_vectors_(state,x,z); 
            }
           
            // Initializes everything
            static void init(
                t& state,
                const X_Vector& x,
                const Y_Vector& y,
                const Z_Vector& z
            ) {
                init_params(state);
                init_vectors(state,x,y,z);
            }

            // Check that we have a valid set of parameters.
            static void check(const Messaging& msg,const t& state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
                InequalityConstrained <Real,XX,ZZ>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Restart();

        public:
            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_nat : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
            struct is_param : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                bool operator () (const std::string& name) const {
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
                const Messaging& msg,
                const Reals& reals,
                const Nats& nats,
                const Params& params,
                const X_Vectors& xs,
                const Y_Vectors& ys,
                const Z_Vectors& zs
            ) {
                peopt::checkLabels <is_real>
                    (msg,reals.first," real name: ");
                peopt::checkLabels <is_nat>
                    (msg,nats.first," natural name: ");
                peopt::checkLabels <is_param>
                    (msg,params.first," paramater name: ");
                peopt::checkLabels <is_x>
                    (msg,xs.first," variable name: ");
                peopt::checkLabels <is_y>
                    (msg,ys.first," equality multiplier name: ");
                peopt::checkLabels <is_z>
                    (msg,zs.first," inequality multiplier name: ");
            }
            
            // Checks whether or not the value used to represent a parameter
            // is valid.  This function returns a string with the error
            // if there is one.  Otherwise, it returns an empty string.
            struct checkParamVal : public std::binary_function
                <std::string,std::string,std::string>
            {
                std::string operator() (
                    std::string label,
                    std::string val
                ) {

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
                typename State::t& state,
                X_Vectors& xs,
                Y_Vectors& ys,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::Restart::stateToVectors(state,xs);
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
                const Messaging& msg,
                typename State::t& state,
                X_Vectors& xs,
                Y_Vectors& ys,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                checkLabels(msg,reals,nats,params,xs,ys,zs);

                // Check the strings used to represent parameters
                checkParams <checkParamVal> (msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::Restart::vectorsToState(state,xs);
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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Functions();

        public:
            // Actual storage of the functions required
            struct t: 
                public EqualityConstrained <Real,XX,YY>::Functions::t,
                public InequalityConstrained <Real,XX,ZZ>::Functions::t
            {
            private:
                // Prevent the use of the copy constructor and the assignment
                // operator.  Since this class holds a number of auto_ptrs
                // to different functions, it is not safe to allow them to
                // be copied.
                t& operator = (const t&);
                t(const t&);

            public:
                
                // Initialize all of the pointers to null
                t() : EqualityConstrained <Real,XX,YY>::Functions::t(), 
                    InequalityConstrained <Real,XX,ZZ>::Functions::t() {}
            };

            // Check that all the functions are defined
            static void check(const Messaging& msg,const t& fns) {
                EqualityConstrained <Real,XX,YY>::Functions::check(msg,fns);
                InequalityConstrained <Real,XX,ZZ>::Functions::check(msg,fns);
            }

            // Initialize any missing functions for just constrained 
            // optimization.
            static void init_(
                const Messaging& msg,
                typename State::t& state,
                t& fns
            ) {
            }

            // Initialize any missing functions 
            static void init(
                const Messaging& msg,
                typename State::t& state,
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
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Printer();

        public:
            // Gets the header for the state information
            static void getStateHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) { 
            }
            // Combines all of the state headers
            static void getStateHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getStateHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer::getStateHeader_
                    (state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer::getStateHeader_
                    (state,out);
            }

            // Gets the state information for output
            static void getState_(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
            }

            // Combines all of the state information
            static void getState(
                const typename Functions::t& fns,
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(fns,state,blank,noiter,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getState_(fns,state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getState_(fns,state,blank,out);
            }
            
            // Get the header for the Krylov iteration
            static void getKrylovHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov headers
            static void getKrylovHeader(
                const typename State::t& state,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylovHeader_(state,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylovHeader_(state,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylovHeader_(state,out);
            }
            
            // Get the information for the Krylov iteration
            static void getKrylov_(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) { }

            // Combines all of the Krylov information
            static void getKrylov(
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer::getKrylov_(state,blank,out);
                EqualityConstrained <Real,XX,YY>::Printer
                    ::getKrylov_(state,blank,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getKrylov_(state,blank,out);
            }
        };
        
        // This contains the different algorithms used for optimization 
        struct Algorithms {
        private:
            // This is a namespace inside of a class.  Do not allow
            // construction.
            Algorithms();

        public:
            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                const Messaging& msg,
                typename Functions::t& fns,
                typename State::t& state
            ){
                // Create an empty state manipulator
                StateManipulator <Constrained <Real,XX,YY,ZZ> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
            
            // Initializes remaining functions then solves an optimization
            // problem
            static void getMin(
                const Messaging& msg,
                const StateManipulator <Constrained <Real,XX,YY,ZZ> >&smanip,
                typename Functions::t& fns,
                typename State::t& state
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
                
                // Minimize the problem
                Unconstrained <Real,XX>::Algorithms
                    ::getMin_(msg,cmanip,fns,state);
            }
        };

    };
}
#endif
