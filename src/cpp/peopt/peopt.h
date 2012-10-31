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
    template <typename Real,template <typename> class XX>
    struct ScalarValuedFunction {
    private:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

    public:
        // <- f(x) 
        virtual Real operator () (const Vector& x) const = 0;

        // g = grad f(x) 
        virtual void grad(const Vector& x,Vector& g) const = 0;

        // H_dx = hess f(x) dx 
        virtual void hessvec(const Vector& x,const Vector& dx,Vector& H_dx)
            const = 0;

        // Allow a derived class to deallocate memory
        virtual ~ScalarValuedFunction() {}
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
    struct AlgorithmClass{
        enum t{
            TrustRegion,            // Trust-Region algorithms
            LineSearch              // Line-search algorithms
        };

        // Converts the algorithm class to a string
        static std::string to_string(t algorithm_class){
            switch(algorithm_class){
            case TrustRegion:
                return "TrustRegion";
            case LineSearch:
                return "LineSearch";
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
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="TrustRegion" ||
                    name=="LineSearch"
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Reasons why we stop the algorithm
    struct StoppingCondition{
        enum t{
            NotConverged,            // Algorithm did not converge
            RelativeGradientSmall,   // Relative gradient was sufficiently small
            RelativeStepSmall,       // Relative change in the step is small
            MaxItersExceeded,        // Maximum number of iterations exceeded
            InteriorPointInstability,// Instability in the interior point method
            External                 // Some external stopping condition 
        };

        // Converts the stopping condition to a string 
        static std::string to_string(t opt_stop){
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
            case External:
                return "External";
            default:
                throw;
            }
        }

        // Converts a string to a stopping condition
        static t from_string(std::string opt_stop){
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
            else if(opt_stop=="External")
                return External;
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
                    name=="External"
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Various operators for both Hessian approximations and preconditioners
    struct Operators{
        enum t{
            Identity,          // Identity approximation
            ScaledIdentity,    // Scaled identity approximation
            BFGS,              // BFGS approximation
            InvBFGS,           // Inverse BFGS approximation
            SR1,               // SR1 approximation
            InvSR1,            // Inverse SR1 approximation
            External           // An external operator provided by the user
        };
        
        // Converts the operator type to a string 
        static std::string to_string(t op){
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
            case External:
                return "External";
            default:
                throw;
            }
        }
        
        // Converts a string to a operator 
        static t from_string(std::string op){
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
            else if(op=="External")
                return External; 
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
                    name=="External" 
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // Different kinds of search directions 
    struct LineSearchDirection{
        enum t{
            SteepestDescent,          // SteepestDescent 
            FletcherReeves,           // Fletcher-Reeves CG
            PolakRibiere,             // Polak-Ribiere CG
            HestenesStiefel,          // HestenesStiefel CG
            BFGS,                     // Limited-memory BFGS 
            NewtonCG                  // Newton-CG
        };
        
        // Converts the line-search direction to a string 
        static std::string to_string(t dir){
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
        static t from_string(std::string dir){
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

    struct LineSearchKind{
        enum t{
            Brents,           // Brent's minimization
            GoldenSection,    // Golden-section search 
            BackTracking,     // BackTracking search 
            TwoPointA,        // Barzilai and Borwein's method A
            TwoPointB         // Barzilai and Borwein's method B
        };
            
        // Converts the line-search kind to a string 
        static std::string to_string(t kind){
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
        static t from_string(std::string kind){
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
    
    struct OptimizationLocation{
        enum t{
            // Occurs before the initial function and gradient evaluation 
            BeforeInitialFuncAndGrad,

            // Occurs after the initial function and gradient evaluation 
            AfterInitialFuncAndGrad,

            // Occurs just before we take the optimization step x+s
            BeforeSaveOld,

            // Occurs just before we take the optimization step x+s
            BeforeStep,

            // Occurs before we calculate our new step.
            BeforeGetStep,

            // Occurs after we take the optimization step u+s, but before
            // we calculate the gradient based on this new step.  In addition,
            // after this point we set the merit value, merit_x, to be
            // merit_xps.
            AfterStepBeforeGradient,

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
        static std::string to_string(t loc){
            switch(loc){
            case BeforeInitialFuncAndGrad:
                return "BeforeInitialFuncAndGrad";
            case AfterInitialFuncAndGrad:
                return "AfterInitialFuncAndGrad";
            case BeforeSaveOld:
                return "BeforeSaveOld";
            case BeforeStep:
                return "BeforeStep";
            case BeforeGetStep:
                return "BeforeGetStep";
            case AfterStepBeforeGradient:
                return "AfterStepBeforeGradient";
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
        static t from_string(std::string loc){
            if(loc=="BeforeInitialFuncAndGrad")
                return BeforeInitialFuncAndGrad;
            else if(loc=="AfterInitialFuncAndGrad")
                return AfterInitialFuncAndGrad;
            else if(loc=="BeforeSaveOld")
                return BeforeSaveOld; 
            else if(loc=="BeforeStep")
                return BeforeStep; 
            else if(loc=="BeforeGetStep")
                return BeforeGetStep; 
            else if(loc=="AfterStepBeforeGradient")
                return AfterStepBeforeGradient; 
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
                if( name=="BeforeInitialFuncAndGrad" ||
                    name=="AfterInitialFuncAndGrad" ||
                    name=="BeforeSaveOld" || 
                    name=="BeforeStep" || 
                    name=="BeforeGetStep" || 
                    name=="AfterStepBeforeGradient" ||
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
    
    struct ProblemClass{
        // Different problem classes that we allow 
        enum t{
            Unconstrained,         // Unconstrained optimization 
            EqualityConstrained,   // Equality constrained optimization 
            InequalityConstrained, // Inequality constrained optimization 
            Constrained            // Fully constrained optimization 
        };

        // Converts the problem class to a string
        static std::string to_string(t problem_class){
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
        static t from_string(std::string problem_class){
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
    
    // Different kinds of interior point methods
    struct InteriorPointMethod{
        enum t{
            PrimalDual,          // Standard primal-dual interior point method 
            PrimalDualLinked,    // A primal dual IPM, but the primal and dual
                                 // variables are kept in lock step.
            LogBarrier,          // Primal log-barrier method 
        };
        
        // Converts the interior point method to a string 
        static std::string to_string(t ipm){
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
        static t from_string(std::string ipm){
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
    struct CentralityStrategy{
        enum t{
            Constant,           // We keep sigma fixed at each iteration.
            PredictorCorrector, // On odd iterations, sigma=1, on even, sigma=0.
        };
        
        // Converts the centrality strategy to a string
        static std::string to_string(t cstrat){
            switch(cstrat){
            case Constant:
                return "Constant";
            case PredictorCorrector:
                return "PredictorCorrector";
            default:
                throw;
            }
        }
        
        // Converts a string to the cstrat
        static t from_string(std::string cstrat){
            if(cstrat=="Constant")
                return Constant; 
            else if(cstrat=="PredictorCorrector")
                return PredictorCorrector; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="Constant" ||
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
        std::ostream& formatReal(std::ostream& out) {
            return out<<std::setprecision(2) << std::scientific << std::setw(10)
                << std::left;
        }
        std::ostream& formatInt(std::ostream& out) {
            return out << std::setw(10) << std::left;
        }
        std::ostream& formatString(std::ostream& out) {
            return out << std::setw(10) << std::left;
        }
        template <typename T>
        std::string atos(T x);
        template <>
        std::string atos <double> (double x){
            std::stringstream ss;
            ss << formatReal << x;
            return ss.str();
        }
        template <>
        std::string atos <unsigned int> (unsigned int x){
            std::stringstream ss;
            ss << formatInt << x;
            return ss.str();
        }
        template <>
        std::string atos <const char*> (const char* x){
            std::stringstream ss;
            ss << formatString << x;
            return ss.str();
        }
        template <>
        std::string atos <KrylovStop::t> (KrylovStop::t x){
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
            default:
                throw;
            }
        }
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
            Real min_rel_err = Real(std::numeric_limits<double>::quiet_NaN());
            for(int i=-2;i<=5;i++){
                Real epsilon=pow(Real(.1),i);
                Real dd=directionalDerivative <> (f,x,dx,epsilon);

                // Calculate the relative error
                Real rel_err=fabs(dd_grad-dd)/(Real(1e-16)+fabs(dd_grad));

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
            Real min_rel_err = Real(std::numeric_limits<double>::quiet_NaN());
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),hess_f_dx,res);

                // Determine the relative error
                Real rel_err=sqrt(X::innr(res,res)) /
                    (Real(1e-16)+sqrt(X::innr(hess_f_dx,hess_f_dx)));

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
            Real min_rel_err = Real(std::numeric_limits<double>::quiet_NaN());
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                Y::axpy(Real(-1.),fp_x_dx,res);

                // Determine the relative error
                Real rel_err=sqrt(Y::innr(res,res)) / 
                    (Real(1e-16)+sqrt(Y::innr(fp_x_dx,fp_x_dx)));
                
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
            Real min_rel_err = Real(std::numeric_limits<double>::quiet_NaN());
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,dy,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),fpps_x_dx_dy,res);

                // Determine the relative error
                Real rel_err=sqrt(X::innr(res,res))
                    / (Real(1e-16)+sqrt(X::innr(fpps_x_dx_dy,fpps_x_dx_dy)));
                
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
        DiagnosticManipulator(
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
            const unsigned int& msg_level=state.msg_level;

            switch(loc){

            // Output the headers for the diagonstic information
            case OptimizationLocation::AfterInitialFuncAndGrad:
                if(msg_level >=1) {
                    // Get the headers 
                    std::list <std::string> out;
                    ProblemClass::Printer::getStateHeader(state,out);
                    if(msg_level >=2)
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
                        ProblemClass::Printer::getState(state,false,true,out);
                    else 
                        ProblemClass::Printer::getState(state,false,false,out);

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
                    ProblemClass::Printer::getState(state,true,false,out);

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
        ConversionManipulator(const StateManipulator <Internal>& smanip_)
            : smanip(smanip_) {}

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
    // where f : X -> R
    template <typename Real,template <typename> class XX> 
    struct Unconstrained {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        typedef std::pair < std::list <std::string>,
                            std::list <Real> > Reals;
        typedef std::pair < std::list <std::string>,
                            std::list <unsigned int> > Nats;
        typedef std::pair < std::list <std::string>,
                            std::list <std::string> > Params; 
        typedef std::pair < std::list <std::string>,
                            std::list <X_Vector> > X_Vectors;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {

            // The actual internal state of the optimization
            struct t {

                // ------------- GENERIC ------------- 

                // Tolerance for the gradient stopping condition
                Real eps_g;

                // Tolerance for the step length stopping criteria
                Real eps_s;

                // Number of control objects to store in a quasi-Newton method
                unsigned int stored_history;

                // Number of failed iterations before we reset the history for
                // quasi-Newton methods
                unsigned int history_reset;

                // Current iteration
                unsigned int iter;

                // Maximum number of optimization iterations
                unsigned int iter_max;

                // Why we've stopped the optimization
                StoppingCondition::t opt_stop;

                // Current number of Krylov iterations taken
                unsigned int krylov_iter;

                // Maximum number of iterations in the Krylov method
                unsigned int krylov_iter_max;

                // Total number of Krylov iterations taken
                unsigned int krylov_iter_total;

                // The maximum number of vectors we orthogonalize against in 
                // the Krylov method.  For something like CG, this is 1.
                unsigned int krylov_orthog_max;

                // Why the Krylov method was last stopped
                KrylovStop::t krylov_stop;

                // Relative error in the Krylov method
                Real krylov_rel_err;

                // Stopping tolerance for the Krylov method
                Real eps_krylov;

                // Algorithm class
                AlgorithmClass::t algorithm_class;

                // Preconditioner
                Operators::t Minv_type;

                // Hessian approximation
                Operators::t H_type;

                // Norm of the gradient
                Real norm_g;

                // Norm of a typical tradient
                Real norm_gtyp;

                // Norm of the trial step
                Real norm_s;

                // Norm of a typical trial step
                Real norm_styp;

                // Optimization variable 
                std::list <X_Vector> x; 
                
                // Gradient 
                std::list <X_Vector> g;
                
                // Trial step 
                std::list <X_Vector> s;
                
                // Old optimization variable 
                std::list <X_Vector> x_old; 
                
                // Old gradient 
                std::list <X_Vector> g_old;
                
                // Old trial step 
                std::list <X_Vector> s_old;

                // Contains the prior iteration information for the
                // quasi-Newton operators
                std::list <X_Vector> oldY;
                std::list <X_Vector> oldS;

                // Current value of the merit function 
                Real merit_x;

                // Merit function at the trial step
                Real merit_xps;

                // Messaging level
                unsigned int msg_level;
                
                // ------------- TRUST-REGION ------------- 

                // Trust region radius
                Real delta;

                // Maximum trust region radius
                Real delta_max;

                // Trust-region parameter for checking whether a step has been
                // accepted
                Real eta1;

                // Trust-region parameter for checking whether a step has been
                // accepted
                Real eta2;

                // Ratio between the predicted and actual reduction
                Real rho;

                // Number of rejected trust-region steps
                unsigned int rejected_trustregion;

                // ------------- LINE-SEARCH ------------- 

                // Line-search step length
                Real alpha;

                // Current number of iterations used in the line-search
                unsigned int linesearch_iter;

                // Maximum number of iterations used in the line-search
                unsigned int linesearch_iter_max;

                // Total number of line-search iterations computed
                unsigned int linesearch_iter_total;

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
                t(const X_Vector& x) {
                    Unconstrained <Real,XX>::State::init_params(*this);
                    Unconstrained <Real,XX>::State::init_vectors(*this,x);
                }
                
                // A trick to allow dynamic casting later
                virtual ~t() {}
            };

            // This sets all of the parameters possible that don't require
            // special memory allocation such as variables.
            static void init_params_(t& state){
                state.eps_g=Real(1e-6);
                state.eps_s=Real(1e-6);
                state.stored_history=0;
                state.history_reset=5;
                state.iter=1;
                state.iter_max=10;
                state.opt_stop=StoppingCondition::NotConverged;
                state.krylov_iter=1;
                state.krylov_iter_max=10;
                state.krylov_iter_total=0;
                state.krylov_orthog_max=1;
                state.krylov_stop=KrylovStop::RelativeErrorSmall;
                state.krylov_rel_err=Real(0.);
                state.eps_krylov=Real(1e-2);
                state.algorithm_class=AlgorithmClass::TrustRegion;
                state.Minv_type=Operators::Identity;
                state.H_type=Operators::Identity;
                state.norm_g=Real(std::numeric_limits<double>::quiet_NaN());
                state.norm_gtyp=Real(std::numeric_limits<double>::quiet_NaN());
                state.norm_s=Real(std::numeric_limits<double>::quiet_NaN());
                state.norm_styp=Real(std::numeric_limits<double>::quiet_NaN());
                state.merit_x=Real(std::numeric_limits<double>::quiet_NaN());
                state.merit_xps=Real(std::numeric_limits<double>::quiet_NaN());
                state.msg_level=1;
                state.delta=Real(100.);
                state.delta_max=Real(100.);
                state.eta1=Real(.1);
                state.eta2=Real(.9);
                state.rho=Real(0.);
                state.rejected_trustregion=0;
                state.alpha=1.;
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
                    X::init(x,state.x.back());
                    X::copy(x,state.x.back());
                state.g.clear();
                    state.g.push_back(X_Vector());
                    X::init(x,state.g.back());
                state.s.clear();
                    state.s.push_back(X_Vector());
                    X::init(x,state.s.back()); 
                state.x_old.clear();
                    state.x_old.push_back(X_Vector());
                    X::init(x,state.x_old.back());
                state.g_old.clear();
                    state.g_old.push_back(X_Vector());
                    X::init(x,state.g_old.back()); 
                state.s_old.clear();
                    state.s_old.push_back(X_Vector());
                    X::init(x,state.s_old.back()); 
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
                if(state.eps_g <= Real(0.)) 
                    ss << "The tolerance for the gradient stopping condition "
                        "must be positive: eps_g = " << state.eps_g;
            
                // Check that the tolerance for the step length stopping
                // condition is positive
                else if(state.eps_s <= Real(0.)) 
                    ss << "The tolerance for the step length stopping "
                        "condition must be positive: eps_s = " << state.eps_s; 
        
                // Check that the current iteration is positive
                else if(state.iter <= 0) 
                    ss << "The current optimization iteration must be "
                        "positive: iter = " << state.iter;

                // Check that the maximum iteration is positive
                else if(state.iter_max <= 0) 
                    ss << "The maximum optimization iteration must be "
                        "positive: iter_max = " << state.iter_max;

                // Check that the current Krylov iteration is positive
                else if(state.krylov_iter <= 0) 
                    ss << "The current Krlov iteration must be "
                        "positive: krylov_iter = " << state.krylov_iter;

                // Check that the maximum Krylov iteration is positive
                else if(state.krylov_iter_max <= 0) 
                    ss << "The maximum Krylov iteration must be "
                        "positive: krylov_iter_max = " << state.krylov_iter_max;

                // Check that the number of vectors we orthogonalize against
                // is at least 1.
                else if(state.krylov_orthog_max <= 0) 
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

                // Check that the norm of the gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_g < Real(0.)
                    || (state.iter!=1 && state.norm_g!=state.norm_g)
                )
                    ss << "The norm of the gradient must be nonnegative: "
                        "norm_g = " << state.norm_g; 

                // Check that the norm of a typical gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_gtyp < Real(0.)
                    || (state.iter!=1 && state.norm_gtyp!=state.norm_gtyp)
                ) 
                    ss << "The norm of a typical gradient must be nonnegative: "
                        "norm_gtyp = " << state.norm_gtyp; 

                // Check that the norm of the trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_s < Real(0.)
                    || (state.iter!=1 && state.norm_s!=state.norm_s)
                ) 
                    ss << "The norm of the trial step must be nonnegative: "
                        "norm_s = " << state.norm_s; 

                // Check that the norm of a typical trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                else if(state.norm_styp < Real(0.)
                    || (state.iter!=1 && state.norm_styp!=state.norm_styp)
                ) 
                    ss << "The norm of a typical trial step must be "
                        "nonnegative: norm_styp = " << state.norm_styp; 

                // Check that the merit functions's value isn't a NaN past
                // iteration 1
                else if(state.iter!=1 && state.merit_x!=state.merit_x)
                    ss<< "The merit function value must be a number: merit_x = "
                        << state.merit_x;

                // Check that the merit function's value at a trial step isn't
                // a NaN past iteration 1
                else if(state.iter!=1 && state.merit_xps!=state.merit_xps) 
                    ss << "The objective value at the trial step must be a "
                        "number: merit_xps = " << state.merit_xps;

                // Check that the trust-region radius is positive
                else if(state.delta<=Real(0.))
                    ss << "The trust-region radius must be positive: delta = "
                        << state.delta; 

                // Check that the maximum trust-region radius is positive
                else if(state.delta_max<=Real(0.))
                    ss << "The maximum trust-region radius must be positive: "
                        "delta_max = " << state.delta_max; 

                // Check that the current trust-region radius is less than
                // or equal to the maximum trust-region radius
                else if(state.delta > state.delta_max)
                    ss << "The trust-region radius must be less than or equal "
                        "to the maximum trust-region radius: delta = "
                        << state.delta << ", delta_max = " << state.delta_max;

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

                // Check that the prediction versus actual reduction is
                // nonnegative 
                else if(state.rho < Real(0.)) 
                    ss << "The predicted versus actual reduction must be "
                        "nonnegative: rho = " << state.rho;

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

            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( name == "eps_g" || 
                        name == "eps_s" || 
                        name == "krylov_rel_err" || 
                        name == "eps_krylov" || 
                        name == "norm_g" || 
                        name == "norm_gtyp" || 
                        name == "norm_s" ||
                        name == "norm_styp" || 
                        name == "merit_x" || 
                        name == "merit_xps" ||
                        name == "delta" || 
                        name == "delta_max" || 
                        name == "eta1" || 
                        name == "eta2" || 
                        name == "rho" || 
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
                    if( name == "algorithm_class" || 
                        name == "opt_stop" || 
                        name == "krylov_stop" ||
                        name == "H_type" || 
                        name == "Minv_type" ||
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
                        name == "g" || 
                        name == "s" || 
                        name == "x_old" || 
                        name == "g_old" || 
                        name == "s_old" || 
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

                    // Check the algorithm class
                    if(label=="algorithm_class") {
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
                    } else if(label=="Minv_type"){
                        if(!Operators::is_valid()(val))
                            ss << base << "preconditioner type: " << val;

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
                xs.first.push_back("g");
                xs.second.splice(xs.second.end(),state.g);
                xs.first.push_back("s");
                xs.second.splice(xs.second.end(),state.s);
                xs.first.push_back("x_old");
                xs.second.splice(xs.second.end(),state.x_old);
                xs.first.push_back("g_old");
                xs.second.splice(xs.second.end(),state.g_old);
                xs.first.push_back("s_old");
                xs.second.splice(xs.second.end(),state.s_old);

                // Write out the quasi-Newton information with sequential names
                {int i=1;
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
                {int i=1;
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
                reals.first.push_back("eps_g");
                reals.second.push_back(state.eps_g);
                reals.first.push_back("eps_s");
                reals.second.push_back(state.eps_s);
                reals.first.push_back("krylov_rel_err");
                reals.second.push_back(state.krylov_rel_err);
                reals.first.push_back("eps_krylov");
                reals.second.push_back(state.eps_krylov);
                reals.first.push_back("norm_g");
                reals.second.push_back(state.norm_g);
                reals.first.push_back("norm_gtyp");
                reals.second.push_back(state.norm_gtyp);
                reals.first.push_back("norm_s");
                reals.second.push_back(state.norm_s);
                reals.first.push_back("norm_styp");
                reals.second.push_back(state.norm_styp);
                reals.first.push_back("merit_x");
                reals.second.push_back(state.merit_x);
                reals.first.push_back("merit_xps");
                reals.second.push_back(state.merit_xps);
                reals.first.push_back("delta");
                reals.second.push_back(state.delta);
                reals.first.push_back("delta_max");
                reals.second.push_back(state.delta_max);
                reals.first.push_back("eta1");
                reals.second.push_back(state.eta1);
                reals.first.push_back("eta2");
                reals.second.push_back(state.eta2);
                reals.first.push_back("rho");
                reals.second.push_back(state.rho);
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
                params.first.push_back("Minv_type");
                params.second.push_back(
                    Operators::to_string(state.Minv_type));
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
                    else if(*name0=="g")
                        state.g.splice(state.g.end(),xs.second,x0);
                    else if(*name0=="s")
                        state.s.splice(state.s.end(),xs.second,x0);
                    else if(*name0=="x_old")
                        state.x_old.splice(state.x_old.end(),xs.second,x0);
                    else if(*name0=="g_old")
                        state.g_old.splice(state.g_old.end(),xs.second,x0);
                    else if(*name0=="s_old")
                        state.s_old.splice(state.s_old.end(),xs.second,x0);
                    else if(*name0=="oldY_")
                        state.oldY.splice(state.oldY.end(),xs.second,x0);
                    else if(*name0=="oldS_")
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
                    if(*name=="eps_g") state.eps_g=*real;
                    else if(*name=="eps_s") state.eps_s=*real;
                    else if(*name=="krylov_rel_err") state.krylov_rel_err=*real;
                    else if(*name=="eps_krylov") state.eps_krylov=*real;
                    else if(*name=="norm_g") state.norm_g=*real;
                    else if(*name=="norm_gtyp") state.norm_gtyp=*real;
                    else if(*name=="norm_s") state.norm_s=*real;
                    else if(*name=="norm_styp") state.norm_styp=*real;
                    else if(*name=="merit_x") state.merit_x=*real;
                    else if(*name=="merit_xps") state.merit_xps=*real;
                    else if(*name=="delta") state.delta=*real;
                    else if(*name=="delta_max") state.delta_max=*real;
                    else if(*name=="eta1") state.eta1=*real;
                    else if(*name=="eta2") state.eta2=*real;
                    else if(*name=="rho") state.rho=*real;
                    else if(*name=="alpha") state.alpha=*real;
                    else if(*name=="eps_ls") state.eps_ls=*real;
                }
            
                // Next, copy in any naturals
                std::list <unsigned int>::iterator nat=nats.second.begin();
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
                    if(*name=="algorithm_class")
                        state.algorithm_class
                            =AlgorithmClass::from_string(*param);
                    else if(*name=="opt_stop")
                        state.opt_stop=StoppingCondition::from_string(*param);
                    else if(*name=="krylov_stop")
                        state.krylov_stop=KrylovStop::from_string(*param);
                    else if(*name=="H_type")
                        state.H_type=Operators::from_string(*param);
                    else if(*name=="Minv_type")
                        state.Minv_type=Operators::from_string(*param);
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
            // Actual storage of the functions required
            struct t{
                // Objective function
                std::auto_ptr <ScalarValuedFunction <Real,XX> > f;

                // Preconditioner for the Hessian of the objective
                std::auto_ptr <Operator <Real,XX,XX> > Minv;

                // Merit function used in globalization
                std::auto_ptr <ScalarValuedFunction <Real,XX> > f_merit;

                // Merit function used for the model of the problem
                std::auto_ptr <ScalarValuedFunction <Real,XX> > model;

                // Symmetric, positive definite operator that changes the
                // shape of the trust-region.
                std::auto_ptr <Operator <Real,XX,XX> > TR_op;
                
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
            // norm(g) / delta_max I.
            class ScaledIdentity : public Operator <Real,XX,XX> {
            private:
                // Norm of the gradient
                const Real& norm_g;

                // Maximum size of the trust-region radius
                const Real& delta_max;
            public:
                ScaledIdentity(const typename State::t& state)
                    : norm_g(state.norm_g), delta_max(state.delta_max) {};

                void operator () (const X_Vector& dx,X_Vector& result) const{
                    X::copy(dx,result);
                    X::scal(norm_g/delta_max,result);
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

                // Messaging device in case the qusi-Newton information is bad
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
                        if(inner_y_s<=0)
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
                    while(1){

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

                // Messaging device in case the qusi-Newton information is bad
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
                    while(1){

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

                // Messaging device in case the qusi-Newton information is bad
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
                        if(inner_y_s<=0)
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
                    int i=0;
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
                            H.reset(new ScaledIdentity (state));
                            break;
                        case Operators::BFGS:
                            H.reset(new BFGS(msg,state));
                            break;
                        case Operators::SR1:
                            H.reset(new SR1(msg,state));
                            break;
                        case Operators::External:
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

                 // g = grad f(x) 
                 void grad(const X_Vector& x,X_Vector& g) const {
                     f->grad(x,g);
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

            // A simple merit function for unconstrained optimization problems.
            // Note, this function does *not* own the memory for the objective.
            // We need to also use the objective in the model, so rather than
            // copying all the memory, we keep the memory for the objective
            // in place and use a reference.
            struct Merit : public peopt::ScalarValuedFunction <Real,XX> {
            private:
                const Messaging& msg;

                // This is a little bit weird since we have a reference to
                // an auto_ptr.  Basically, the Functions::t structure owns
                // the memory for all of the functions, which are stored
                // as auto_ptrs.  We want this function to have a reference
                // to that objective.  Now, during the initialization, the
                // objective may and probably be modified for a variety of
                // reasons.  For example, we may replace the function with
                // another that has a modified Hessian.  When this happens,
                // the underlying function is replaced, so if we just have 
                // a reference to the actual function and not the pointer,
                // we could get into trouble.
                const std::auto_ptr <peopt::ScalarValuedFunction <Real,XX> >& f;

                // Make sure we call the constructor below
                Merit () {}
            public:
                Merit(
                    const Messaging& msg_,
                    const typename Functions::t& fns
                ) : msg(msg_), f(fns.f) {}

                // <- f(x) 
                Real operator () (const X_Vector& x) const {
                    return (*f)(x); 
                }

                // Throw an error 
                void grad(const X_Vector& x,X_Vector& g) const {
                    msg.error("The gradient of the merit function is "
                        "undefined.");
                }

                // Throw an error 
                void hessvec(
                    const X_Vector& x,
                    const X_Vector& dx,
                    X_Vector& H_dx
                ) const {
                    msg.error("The Hessian of the merit function is "
                        "undefined.");
                }
            };
           
            // A model for an unconstrained optimization problem.  In our case,
            // we use a quadratic model based on a Taylor series.  
            struct Model: public peopt::ScalarValuedFunction <Real,XX> {
            private:
                const Messaging& msg;
                const peopt::ScalarValuedFunction <Real,XX>& f;
                const Real& merit_x;
                const std::list <X_Vector>& gg;
                const std::list <X_Vector>& xx;
            public:
                Model(
                    const Messaging& msg_,
                    const typename Functions::t& fns,
                    const typename State::t& state
                ) : msg(msg_), f(*(fns.f)), merit_x(state.merit_x), gg(state.g),
                    xx(state.x)
                { }

                // <- m(s) 
                Real operator () (const X_Vector& s) const {
                    // Get the current iterate and the gradient
                    const X_Vector& g=*(gg.begin());
                    const X_Vector& x=*(xx.begin());

                    // Determine H(x)s
                    X_Vector Hx_s; X::init(x,Hx_s);
                    f.hessvec(x,s,Hx_s);

                    // Calculate the model
                    return merit_x+X::innr(g,s)+Real(.5)*X::innr(Hx_s,s);
                }

                // Throw an error 
                void grad(const X_Vector& x,X_Vector& g) const {
                    msg.error("The gradient of the merit function is "
                        "undefined.");
                }

                // Throw an error 
                void hessvec(
                    const X_Vector& x,
                    const X_Vector& dx,
                    X_Vector& H_dx
                ) const {
                    msg.error("The Hessian of the merit function is "
                        "undefined.");
                }
            };

            // Check that all the functions are defined
            static void check(const Messaging& msg,const t& fns) {
                // Check that objective function exists 
                if(fns.f.get()==NULL)
                    msg.error("Missing an objective function definition.");
                
                // Check that a preconditioner exists 
                if(fns.Minv.get()==NULL)
                    msg.error("Missing a preconditioner definition.");
            }

            // Initialize any missing functions for just unconstrained
            // optimization.
            static void init_(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
            ) {
                // Determine the preconditioner
                switch(state.Minv_type){
                    case Operators::Identity:
                        fns.Minv.reset(new Identity());
                        break;
                    case Operators::InvBFGS:
                        fns.Minv.reset(new InvBFGS(msg,state));
                        break;
                    case Operators::InvSR1:
                        fns.Minv.reset(new InvSR1(msg,state));
                        break;
                    case Operators::External:
                        if(fns.Minv.get()==NULL)
                            msg.error("An externally defined preconditioner "
                                "must be provided explicitely.");
                        break;
                    default:
                        msg.error("Not a valid Hessian approximation.");
                        break;
                }

                // If a trust-region operator has not been provided, use the
                // identity.
                if(fns.TR_op.get()==NULL)
                    fns.TR_op.reset(new Identity());

                // Check that all functions are defined (namely, the 
                // objective).
                check(msg,fns);

                // Modify the objective function if necessary
                fns.f.reset(new HessianAdjustedFunction(msg,state,fns));

                // Create the merit function
                fns.f_merit.reset(new Merit(msg,fns));

                // Create the model
                fns.model.reset(new Model(msg,fns,state));
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
                out.push_back(atos <> ("Merit(x)"));
                out.push_back(atos <> ("|| g ||"));
                out.push_back(atos <> ("|| s ||"));

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    out.push_back(atos <> ("Kry Iter"));
                    out.push_back(atos <> ("Kry Err"));
                    out.push_back(atos <> ("Kry Why"));
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    out.push_back(atos <> ("LS Iter"));
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
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {

                // Create some shortcuts
                const unsigned int& iter=state.iter;
                const Real& merit_x=state.merit_x;
                const Real& norm_g=state.norm_g;
                const Real& norm_s=state.norm_s;
                const unsigned int& krylov_iter=state.krylov_iter;
                const Real& krylov_rel_err=state.krylov_rel_err;
                const KrylovStop::t& krylov_stop=state.krylov_stop;
                const unsigned int& linesearch_iter=state.linesearch_iter;
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;
                const unsigned int& rejected_trustregion
                    =state.rejected_trustregion;

                // Figure out if we're at the absolute beginning of the
                // optimization.  We have to be a little saavy about this
                // since we could be on the first iteration, but in the
                // middle of a line-search or trust-region method and
                // still want to output things
                bool opt_begin = (iter==1) &&
                    ((algorithm_class == AlgorithmClass::LineSearch && 
                        linesearch_iter==0) ||
                    (algorithm_class == AlgorithmClass::TrustRegion && 
                        rejected_trustregion == 0));


                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Basic information
                if(!noiter)
                    out.push_back(atos <> (iter));
                else
                    out.push_back(atos <> ("*"));
                out.push_back(atos <> (merit_x));
                out.push_back(atos <> (norm_g));
                if(!opt_begin)
                    out.push_back(atos <> (norm_s));
                else
                    out.push_back("          ");

                // In case we're using a Krylov method
                if(    algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    if(!opt_begin) {
                        out.push_back(atos <> (krylov_iter));
                        out.push_back(atos <> (krylov_rel_err));
                        out.push_back(atos <> (krylov_stop));
                    } else 
                        for(int i=0;i<3;i++) out.push_back("          ");
                }

                // In case we're using a line-search method
                if(algorithm_class==AlgorithmClass::LineSearch) {
                    if(!opt_begin)
                        out.push_back(atos <> (linesearch_iter));
                    else 
                        out.push_back("          ");
                }

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)="          ";
            }

            // Combines all of the state information
            static void getState(
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(state,blank,noiter,out);
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
                const unsigned int& krylov_iter=state.krylov_iter;
                const unsigned int& krylov_iter_total=state.krylov_iter_total;
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

            // Checks a set of stopping conditions
            static StoppingCondition::t checkStop(
                const typename State::t& state
            ){
                // Create some shortcuts
                const Real& norm_g=state.norm_g;
                const Real& norm_gtyp=state.norm_gtyp;
                const Real& norm_s=state.norm_s;
                const Real& norm_styp=state.norm_styp;
                const int& iter=state.iter;
                const int& iter_max=state.iter_max;
                const Real& eps_g=state.eps_g;
                const Real& eps_s=state.eps_s;

                // Check if we've exceeded the number of iterations
                if(iter>=iter_max)
                    return StoppingCondition::MaxItersExceeded;

                // Check whether the change in the step length has become too
                // small relative to some typical step
                if(norm_s < eps_s*norm_styp)
                    return StoppingCondition::RelativeStepSmall;
                
                // Check whether the norm is small relative to some typical
                // gradient
                if(norm_g < eps_g*norm_gtyp)
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

                // Store a reference to the base of the Hessian-vector product
                const X_Vector& x;

            public:
                // Take in the objective and the base point during construction 
                HessianOperator(
                    const ScalarValuedFunction <Real,XX>& f_,
                    const X_Vector& x_) : f(f_), x(x_) {}

                // Basic application
                void operator () (const X_Vector& dx,X_Vector &y) const {
                    f.hessvec(x,dx,y);
                }
            };
        
            // Checks whether we accept or reject a step
            static bool checkStep(
                const typename Functions::t& fns,
                typename State::t& state
            ){
                // Create shortcuts to some elements in the state
                const X_Vector& s=*(state.s.begin());
                const X_Vector& x=*(state.x.begin());
                const Real& eta1=state.eta1;
                const Real& eta2=state.eta2;
                const Real& delta_max=state.delta_max;
                const Real& merit_x=state.merit_x;
                const Real& norm_s=state.norm_s;
                Real& delta=state.delta;
                Real& rho=state.rho;
                Real& merit_xps=state.merit_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f_merit=*(fns.f_merit);
                const ScalarValuedFunction <Real,XX>& model=*(fns.model);

                // Allocate memory for temporaries that we need
                X_Vector xps; X::init(x,xps);

                // Determine x+s 
                X::copy(s,xps);
                X::axpy(Real(1.),x,xps);

                // Determine the merit function evaluated at x+s
                merit_xps=f_merit(xps);
                
                // Determine the model function at s, m(s)
                Real model_s=model(s);

                // Add a safety check in case we don't actually minimize the TR
                // subproblem correctly. This could happen for a variety of
                // reasons.  Most notably, if we do not correctly calculate the
                // Hessian approximation, we could have a nonsymmetric 
                // approximation.  In that case, truncated-CG will exit, but 
                // has an undefined result.  In the case that the actual 
                // reduction also increases, rho could have an extraneous 
                // positive value.  Hence, we require an extra check.
                if(model_s > merit_x){
                    delta = norm_s/Real(2.);
                    rho = Real(std::numeric_limits<double>::quiet_NaN()); 
                    return false;
                }

                // Determine the ratio of reductions
                rho = (merit_x - merit_xps) / (merit_x - model_s);

                // Update the trust region radius and return whether or not we
                // accept the step
                if(rho >= eta2){
                    // Only increase the size of the trust region if we were
                    // close to the boundary
                    if(fabs(norm_s-delta) < Real(1e-4)*delta)
                        delta = std::min(delta*Real(2.),delta_max);
                    return true;
                } else if(rho >= eta1 && rho < eta2)
                    return true;
                else {
                    delta = norm_s/Real(2.);
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
                const Real& eps_s=state.eps_s;
                const Real& eps_krylov=state.eps_krylov;
                const unsigned int& krylov_iter_max=state.krylov_iter_max;
                const unsigned int& krylov_orthog_max=state.krylov_orthog_max;
                const Real& delta=state.delta;
                const X_Vector& x=state.x.back();
                const X_Vector& g=state.g.back();
                const Real& norm_g=state.norm_g;
                const Real& norm_styp=state.norm_styp;
                unsigned int& rejected_trustregion=state.rejected_trustregion;
                X_Vector& s=state.s.back();
                Real& norm_s=state.norm_s;
                unsigned int& krylov_iter=state.krylov_iter;
                Real& krylov_rel_err=state.krylov_rel_err;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                std::list <X_Vector>& oldY=state.oldY; 
                std::list <X_Vector>& oldS=state.oldS; 
                unsigned int& history_reset=state.history_reset;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const Operator <Real,XX,XX>& Minv=*(fns.Minv);
                const Operator <Real,XX,XX>& TR_op=*(fns.TR_op);

                // Allocate some memory for the scaled trial step
                X_Vector x_tmp1; X::init(x,x_tmp1);

                // Continue to look for a step until one comes back as valid
                for(rejected_trustregion=0;
                    true; 
                ) {
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeGetStep);

                    // Use truncated-cd to find a new trial step
                    HessianOperator H(f,x);
                    X_Vector s_cp; X::init(x,s_cp);
                    X_Vector minus_g; X::init(x,minus_g);
                    X::copy(g,minus_g); X::scal(Real(-1.),minus_g);
                    X::zero(s);
                    truncated_pcd(
                        H,
                        minus_g,
                        Minv,
                        TR_op,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        delta,
                        s,
                        s_cp,
                        krylov_rel_err,
                        krylov_iter,
                        krylov_stop);
                    krylov_rel_err = krylov_rel_err / (Real(1e-16)+norm_g);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::BeforeActualVersusPredicted);

                    // Save the length of the scaled trial step
                    TR_op(s,x_tmp1);
                    norm_s=sqrt(X::innr(x_tmp1,x_tmp1));

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
                    if(norm_s < eps_s*norm_styp) {
                        X::scal(Real(0.),s);
                        norm_s=Real(0.);
                        break;
                    }
                } 
            }
        
            // Steepest descent search direction
            static void SteepestDescent(
                typename State::t& state
            ) {
                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                X_Vector& s=*(state.s.begin());

                // We take the steepest descent direction
                X::copy(g,s);
                X::scal(Real(-1.),s);
            }
    
            // Nonlinear Conjugate Gradient
            static void NonlinearCG(
                const NonlinearCGDirections::t dir,
                typename State::t& state
            ) {
            
                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                const X_Vector& s_old=*(state.s_old.begin());
                const int& iter=state.iter;
                X_Vector& s=*(state.s.begin());

                // If we're on the first iterations, we take the steepest
                // descent direction
                if(iter==1) SteepestDescent(state);

                // On subsequent iterations, we take the specified direction
                else {
                    // Find the momentum parameter
                    double beta=0;
                    switch(dir) {
                    case NonlinearCGDirections::FletcherReeves:
                        beta=FletcherReeves(state);
                        break;
                    case NonlinearCGDirections::PolakRibiere:
                        beta=PolakRibiere(state);
                        break;
                    case NonlinearCGDirections::HestenesStiefel:
                        beta=HestenesStiefel(state);
                        break;
                    }

                    // Find -g+beta*s_old
                    X::copy(g,s);
                    X::scal(Real(-1.),s);
                    X::axpy(beta,s_old,s);
                }
            }

            // Fletcher-Reeves CG search direction
            static Real FletcherReeves(
                typename State::t& state
            ) {

                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                const X_Vector& g_old=*(state.g_old.begin());

                // Return the momentum parameter
                return X::innr(g,g)/X::innr(g_old,g_old);
            }
        
            // Polak-Ribiere CG search direction
            static Real PolakRibiere(
                typename State::t& state
            ) {

                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                const X_Vector& g_old=*(state.g_old.begin());
                    
                // Return the momentum parameter
                return (X::innr(g,g)-X::innr(g,g_old))
                    /X::innr(g_old,g_old);
            }
            
            // Hestenes-Stiefel search direction
            static Real HestenesStiefel(
                typename State::t& state
            ) {

                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                const X_Vector& g_old=*(state.g_old.begin());
                const X_Vector& s_old=*(state.s_old.begin());
                    
                // Return the momentum parameter
                return (X::innr(g,g)-X::innr(g,g_old))
                    /(X::innr(g,s_old)-X::innr(g_old,s_old));
            }

            // BFGS search direction
            static void BFGS(
                const Messaging& msg,
                typename State::t& state
            ) {
                
                // Create some shortcuts 
                const X_Vector& g=*(state.g.begin());
                X_Vector& s=*(state.s.begin());

                // Create the inverse BFGS operator
                typename Functions::InvBFGS Hinv(msg,state); 

                // Apply the inverse BFGS operator to the gradient
                Hinv(g,s);

                // Negate the result
                X::scal(Real(-1.),s);
            }

            // Compute a Golden-Section search between eps and 2*alpha where
            // alpha is the last line search parameter.
            static void goldenSection(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const X_Vector& x=*(state.x.begin());
                const unsigned int& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                X_Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& merit_xps=state.merit_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f_merit=*(fns.f_merit);

                // Create one work element that holds x+mu s or x+lambda s
                X_Vector x_p_s; X::init(x,x_p_s);

                // Find 1 over the golden ratio
                Real beta=Real(2./(1.+sqrt(5.)));

                // Find a bracket for the linesearch such that a < b
                Real a=Real(0.);
                Real b=Real(2.)*alpha;

                // Find two new points between a and b, mu and lambda,
                // such that lambda < mu
                double lambda=a+(1.-beta)*(b-a);
                double mu=a+beta*(b-a);

                // Find the merit value at mu and labmda 

                // mu 
                X::copy(x,x_p_s);
                X::axpy(mu,s,x_p_s);
                Real merit_mu=f_merit(x_p_s);

                // lambda
                X::copy(x,x_p_s);
                X::axpy(lambda,s,x_p_s);
                Real merit_lambda=f_merit(x_p_s);

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

                        X::copy(x,x_p_s);
                        X::axpy(mu,s,x_p_s);
                        merit_mu=f_merit(x_p_s);

                    // Otherwise, the objective is greater on the right, so
                    // bracket on the left
                    } else {
                        b=mu;
                        mu=lambda;
                        merit_mu=merit_lambda;
                        lambda=a+(1-beta)*(b-a);
                
                        X::copy(x,x_p_s);
                        X::axpy(lambda,s,x_p_s);
                        merit_lambda=f_merit(x_p_s);
                    }
                }

                // Keep track of the total number of iterations
                iter_total += iter;

                // Once we're finished narrowing in on a solution, take our best
                // guess for the line search parameter
                alpha=merit_lambda < merit_mu ? lambda : mu;

                // Save the objective value at this step
                merit_xps=merit_lambda < merit_mu ? merit_lambda : merit_mu;
            }

            // Find the line search parameter based on the 2-point approximation
            // from Barzilai and Borwein
            static void twoPoint(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const X_Vector& x=*(state.x.begin());
                const X_Vector& g=*(state.g.begin());
                const X_Vector& x_old=*(state.x_old.begin());
                const X_Vector& g_old=*(state.g_old.begin());
                const LineSearchKind::t& kind=state.kind;
                Real& alpha=state.alpha;
                X_Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& merit_xps=state.merit_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f_merit=*(fns.f_merit);

                // Create elements for delta_x and delta_g as well as one work
                // element for storing x+alpha s
                X_Vector delta_x; X::init(x,delta_x);
                X_Vector delta_g; X::init(x,delta_g);
                X_Vector x_p_s; X::init(x,x_p_s);

                // Find delta_x
                X::copy(x,delta_x);
                X::axpy(Real(-1.),x_old,delta_x);

                // Find delta_g
                X::copy(g,delta_g);
                X::axpy(Real(-1.),g_old,delta_g);

                // Find alpha
                if(kind==LineSearchKind::TwoPointA)
                    alpha=X::innr(delta_x,delta_g)/X::innr(delta_g,delta_g);
                else if(kind==LineSearchKind::TwoPointB)
                    alpha=X::innr(delta_x,delta_x)/X::innr(delta_x,delta_g);

                // Save the merit value at this step
                X::copy(x,x_p_s);
                X::axpy(alpha,s,x_p_s);
                merit_xps=f_merit(x_p_s);

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
                const X_Vector& x=*(state.x.begin());
                const unsigned int& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                X_Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& merit_xps=state.merit_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f_merit=*(fns.f_merit);

                // Create one work element for holding x+alpha s
                X_Vector x_p_s; X::init(x,x_p_s);

                // Store the best merit value and alpha that we used to find it.
                // Our initial guess will be at alpha*2.
                Real alpha_best=Real(2.)*alpha;
                X::copy(x,x_p_s);
                X::axpy(alpha_best,s,x_p_s);
                Real merit_best=f_merit(x_p_s);

                // Evaluate the merit iter_max times at a distance of
                // 2*alpha, alpha, alpha/2, ....  Then, pick the best one.
                // Note, we start iter at 1 since we've already done one
                // iteration above.
                Real alpha0=alpha;
                for(iter=1;iter<iter_max;iter++){
                    // Evaluate f_merit(x+alpha*s)
                    X::copy(x,x_p_s);
                    X::axpy(alpha0,s,x_p_s);
                    Real merit=f_merit(x_p_s);

                    // If this is better than our best guess so far, save it
                    if(merit<merit_best){
                        merit_best=merit;
                        alpha_best=alpha0;
                    }

                    // Reduce the size of alpha
                    alpha0 /= Real(2.);
                }

                // Save the best merit value and alpha found
                alpha=alpha_best;
                merit_xps=merit_best;

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
                const X_Vector& x=state.x.back();
                const X_Vector& g=state.g.back();
                const LineSearchDirection::t& dir=state.dir;
                const LineSearchKind::t& kind=state.kind;
                const int& iter=state.iter;
                const int& linesearch_iter_max=state.linesearch_iter_max;
                const Real& merit_x=state.merit_x;
                const Real& eps_s=state.eps_s;
                const Real& norm_styp=state.norm_styp;
                const Real& eps_krylov=state.eps_krylov;
                const unsigned int& krylov_iter_max=state.krylov_iter_max;
                const unsigned int& krylov_orthog_max=state.krylov_orthog_max;
                const Real& norm_g=state.norm_g;
                X_Vector& s=state.s.back();
                Real& merit_xps=state.merit_xps;
                Real& norm_s=state.norm_s;
                Real& alpha=state.alpha;
                Real& krylov_rel_err=state.krylov_rel_err;
                unsigned int& krylov_iter=state.krylov_iter;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const Operator <Real,XX,XX>& Minv=*(fns.Minv);
                const Operator <Real,XX,XX>& TR_op=*(fns.TR_op);
                    
                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeGetStep);

                // Find the line-search direction
                switch(dir){
                case LineSearchDirection::SteepestDescent:
                    SteepestDescent(state);
                    break;
                case LineSearchDirection::FletcherReeves:
                    NonlinearCG(NonlinearCGDirections::FletcherReeves,state);
                    break;
                case LineSearchDirection::PolakRibiere:
                    NonlinearCG(NonlinearCGDirections::PolakRibiere,state);
                    break;
                case LineSearchDirection::HestenesStiefel:
                    NonlinearCG(NonlinearCGDirections::HestenesStiefel,state);
                    break;
                case LineSearchDirection::BFGS:
                    BFGS(msg,state);
                    break;
                case LineSearchDirection::NewtonCG: {
                    HessianOperator H(f,x);
                    X_Vector s_cp; X::init(x,s_cp);
                    X_Vector minus_g; X::init(x,minus_g);
                    X::copy(g,minus_g); X::scal(Real(-1.),minus_g);
                    X::zero(s);
                    truncated_pcd(
                        H,
                        minus_g,
                        Minv,
                        TR_op,
                        eps_krylov,
                        krylov_iter_max,
                        krylov_orthog_max,
                        std::numeric_limits <Real>::infinity(),
                        s,
                        s_cp,
                        krylov_rel_err,
                        krylov_iter,
                        krylov_stop);
                    krylov_rel_err = krylov_rel_err / (Real(1e-16)+norm_g);
                    break;
                }}
                    
                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::BeforeLineSearch);

                // Do a line-search in the specified direction
                switch(kind){
                case LineSearchKind::GoldenSection:
                    // Continue doing a line-search until we get a reduction
                    // in the merit value.
                    do {
                        // Conduct the golden section search
                        goldenSection(fns,state);

                        // If we have no reduction in the merit, print
                        // some diagnostic information.
                        if(merit_xps > merit_x || merit_xps!=merit_xps) {

                            // Determine the size of the step
                            norm_s=alpha*sqrt(X::innr(s,s));

                            // Check if the step becomes so small that we're not
                            // making progress.  In this case, take a zero step 
                            // and allow the stopping conditions to exit
                            if(norm_s < eps_s*norm_styp) {
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
                    } while(merit_x < merit_xps || merit_xps!=merit_xps);
                    break;
                case LineSearchKind::BackTracking:
                    // Continue doing a line-search until we get a reduction
                    // in the merit value.
                    do {
                        // Conduct the backtracking search
                        backTracking(fns,state);

                        // If we have no reduction in the merit, print
                        // some diagnostic information.
                        if(merit_xps > merit_x || merit_xps!=merit_xps) {
                            // Determine the size of the step
                            norm_s=alpha*sqrt(X::innr(s,s));

                            // Check if the step becomes so small that we're not
                            // making progress.  In this case, take a zero step 
                            // and allow the stopping conditions to exit
                            if(norm_s < eps_s*norm_styp) {
                                alpha=0.;
                                break;
                            }

                            // Manipulate the state if required
                            smanip(fns,state,
                                OptimizationLocation::AfterRejectedLineSearch);

                            // We set alpha to be four times less than the
                            // minimimum alpha we searched before.  We do this
                            // since the line-search always looks twice alpha
                            // out in the beginning of the search.
                            alpha = alpha/pow(Real(2.),linesearch_iter_max+1);
                        }

                    // If we don't decrease the merit, try again 
                    } while(merit_x < merit_xps || merit_xps!=merit_xps);
                    break;
                case LineSearchKind::TwoPointA:
                case LineSearchKind::TwoPointB:
                    if(iter>1) twoPoint(fns,state);
                    else goldenSection(fns,state);
                    break;
                case LineSearchKind::Brents:
                    msg.error(
                        "Brent's linesearch is not currently implemented.");
                    break;
                }
            
                // Scale the line-search direction by the line search parameter 
                X::scal(alpha,s);

                // Save the length of the trial step
                norm_s=sqrt(X::innr(s,s));
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
                }
            }

            // Updates the quasi-Newton information
            static void updateQuasi(
                typename State::t& state
            ){
                // Exit immediately if we're not using a quasi-Newton method
                if(state.stored_history==0) return;

                // Create some shortcuts
                const X_Vector& x=*(state.x.begin());
                const X_Vector& g=*(state.g.begin());
                const X_Vector& x_old=*(state.x_old.begin());
                const X_Vector& g_old=*(state.g_old.begin());
                const Operators::t& Minv_type=state.Minv_type;
                const Operators::t& H_type=state.H_type;
                const LineSearchDirection::t& dir=state.dir;
                std::list <X_Vector>& oldY=state.oldY;
                std::list <X_Vector>& oldS=state.oldS;
               
                // Allocate some temp storage for y and s
                X_Vector s; X::init(x,s);
                X_Vector y; X::init(x,y);

                // Find s = u-u_old
                X::copy(x,s);
                X::axpy(Real(-1.),x_old,s);

                // Find y = g - g_old
                X::copy(g,y);
                X::axpy(Real(-1.),g_old,y);

                // If we're using BFGS, check that <y,s> > 0
                if((Minv_type==Operators::InvBFGS || H_type==Operators::BFGS
                    || dir==LineSearchDirection::BFGS)
                    && X::innr(y,s) <= 0)
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
                X_Vector& x=*(state.x.begin());
                X_Vector& g=*(state.g.begin());
                X_Vector& s=*(state.s.begin());
                X_Vector& x_old=*(state.x_old.begin());
                X_Vector& g_old=*(state.g_old.begin());
                X_Vector& s_old=*(state.s_old.begin());
                Real& merit_x=state.merit_x;
                Real& merit_xps=state.merit_xps;
                Real& norm_s=state.norm_s;
                Real& norm_g=state.norm_g;
                Real& norm_gtyp=state.norm_gtyp;
                Real& norm_styp=state.norm_styp;
                unsigned int& iter=state.iter;
                StoppingCondition::t& opt_stop=state.opt_stop;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const ScalarValuedFunction <Real,XX>& f_merit=*(fns.f_merit);

                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation
                    ::BeforeInitialFuncAndGrad);

                // Evaluate the merit function and gradient if we've not
                // done so already
                if(merit_x != merit_x) {
                    f.grad(x,g);
                    merit_x=f_merit(x);
                    norm_g=sqrt(X::innr(g,g));
                    norm_gtyp=norm_g;

                    // This one is a little funny.  Sometimes, we run into
                    // trouble trying to calculate the initial step.  In this
                    // case, we don't know when to exit due to the relative
                    // step size being small since we've not calculated the
                    // typical step yet.  This safeguards this case.  Basically,
                    // it initially sets the norm of a typical step to be the
                    // norm of the gradient, which is akin to taking a
                    // steepest descent step without globalization.
                    norm_styp=norm_g;
                }

                // Manipulate the state if required
                smanip(fns,state,OptimizationLocation::AfterInitialFuncAndGrad);

                // Primary optimization loop
                do{
                    // Get a new optimization iterate.  
                    getStep(msg,smanip,fns,state);

                    // If we've not calculated it already, save the size of
                    // the step
                    if(norm_styp!=norm_styp) norm_styp=norm_s;
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeSaveOld);

                    // Save the old variable, gradient, and trial step.  This
                    // is useful for both CG and quasi-Newton methods.
                    X::copy(x,x_old);
                    X::copy(g,g_old);
                    X::copy(s,s_old);

                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeStep);

                    // Move to the new iterate
                    X::axpy(Real(1.),s,x);

                    // Manipulate the state if required
                    smanip(fns,state,
                        OptimizationLocation::AfterStepBeforeGradient);

                    // Find the new merit value and gradient
                    merit_x=merit_xps;
                    f.grad(x,g);
                    norm_g=sqrt(X::innr(g,g));
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::BeforeQuasi);

                    // Update the quasi-Newton information
                    updateQuasi(state);
                    
                    // Manipulate the state if required
                    smanip(fns,state,OptimizationLocation::AfterQuasi);

                    // Increase the iteration
                    iter++;
                    
                    // Check the stopping condition
                    opt_stop=checkStop(state);

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
        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector;
        
        typedef std::pair < std::list <std::string>,
                            std::list <Real> > Reals;
        typedef std::pair < std::list <std::string>,
                            std::list <unsigned int> > Nats;
        typedef std::pair < std::list <std::string>,
                            std::list <std::string> > Params; 
        typedef std::pair < std::list <std::string>,
                            std::list <X_Vector> > X_Vectors;
        typedef std::pair < std::list <std::string>,
                            std::list <Y_Vector> > Y_Vectors;

        // Routines that manipulate the internal state of the optimization 
        // algorithm.
        struct State {

            // The actual internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {

                // The Lagrange multiplier (dual variable) for the equality
                // constraints
                std::list <Y_Vector> y;

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
            static void init_params_(t& state) { }
            static void init_params(t& state) {
                Unconstrained <Real,XX>::State::init_params_(state); 
                EqualityConstrained <Real,XX,YY>::State::init_params_(state);
            }

            // This initializes all the variables required for equality
            // constrained optimization.  
            static void init_vectors_(t& state, const Y_Vector& y) {
                state.y.clear();
                    state.y.push_back(Y_Vector());
                    Y::init(y,state.y.back());
                    Y::copy(y,state.y.back());
            }
            static void init_vectors(
                t& state,
                const X_Vector& x,
                const Y_Vector& y
            ) {
                Unconstrained <Real,XX>::State::init_vectors_(state,x); 
                EqualityConstrained <Real,XX,YY>::State::init_vectors_(state,y);
            }
           
            // Initializes everything
            static void init(t& state, const X_Vector& x, const Y_Vector& y) {
                init_params(state);
                init_vectors(state,x,y);
            }

            // Check that we have a valid set of parameters.  
            static void check_(const Messaging& msg,const t& state) { }
            static void check(const Messaging& msg,const t& state) {
                Unconstrained <Real,XX>::State::check_(msg,state);
                EqualityConstrained <Real,XX,YY>::State::check_(msg,state);
            }
        };
        
        // Utilities for restarting the optimization
        struct Restart {

            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( typename Unconstrained <Real,XX>::Restart
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
                    if( typename Unconstrained <Real,XX>::Restart
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
                    if( name == "y") 
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
                    }
                    return ss.str();
                }
            };
            
            // Copy out all equality multipliers 
            static void stateToVectors(
                typename State::t& state, 
                Y_Vectors& ys
            ) {
                ys.first.push_back("y");
                ys.second.splice(ys.second.end(),state.y);
            }

            // Copy out all the scalar information
            static void stateToScalars(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }
            
            // Copy in all equality multipliers 
            static void vectorsToState(
                typename State::t& state,
                Y_Vectors& ys
            ) {
                for(std::list <std::string>::iterator name=ys.first.begin();
                    name!=ys.first.end();
                    name++
                ) {
                    // Since we're using a splice operation, we slowly empty
                    // the multiplier list.  Hence, we always take the first
                    // element.
                    typename std::list <Y_Vector>::iterator yy
                        =ys.second.begin();

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name=="y") state.y.splice(state.y.end(),ys.second,yy);
                }
            }
            
            // Copy in all the scalar information
            static void scalarsToState(
                typename State::t& state,
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }

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
                    ::Restart::stateToVectors(state,ys);
            
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
                    ::Restart::vectorsToState(state,ys);
                
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

            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Equality constraints 
                std::auto_ptr <VectorValuedFunction <Real,XX,YY> > g;
            };

            // Check that all the functions are defined
            static void check(const Messaging& msg,const t& fns) {

                // Check the unconstrained pieces
                Unconstrained <Real,XX>::Functions::check(msg,fns);
                
                // Check that the equality constraints exist 
                if(fns.g.get()==NULL)
                    msg.error("Missing the equality constraint definition.");
            }

            // Initialize any missing functions for just equality constrained 
            // optimization.
            static void init_(
                const Messaging& msg,
                const typename State::t& state,
                t& fns
            ) {
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
        
        // Create some shortcuts for some type names
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        typedef ZZ <Real> Z;
        typedef typename Z::Vector Z_Vector;
        
        typedef std::pair < std::list <std::string>,
                            std::list <Real> > Reals;
        typedef std::pair < std::list <std::string>,
                            std::list <unsigned int> > Nats;
        typedef std::pair < std::list <std::string>,
                            std::list <std::string> > Params; 
        typedef std::pair < std::list <std::string>,
                            std::list <X_Vector> > X_Vectors;
        typedef std::pair < std::list <std::string>,
                            std::list <Z_Vector> > Z_Vectors;

        // Functions that manipulate the internal state of the optimization 
        // algorithm.
        struct State {

            // The actual internal state of the optimization
            struct t: public virtual Unconstrained <Real,XX>::State::t {
                // Original gradient (possibly of a Lagrangian) 
                std::list <X_Vector> g_orig;

                // Gradient used in the Schur complement system. 
                std::list <X_Vector> g_schur;

                // Gradient of the inequality Lagrangian.
                std::list <X_Vector> g_lag;

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
                state.mu = Real(std::numeric_limits<double>::quiet_NaN());
                state.mu_typ = Real(std::numeric_limits<double>::quiet_NaN());
                state.eps_mu= Real(1e-6);
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
                // Allocate memory for g_orig
                state.g_orig.clear();
                    state.g_orig.push_back(X_Vector());
                    X::init(x,state.g_orig.back());

                // Allocate memory for g_schur 
                state.g_schur.clear();
                    state.g_schur.push_back(X_Vector());
                    X::init(x,state.g_schur.back());

                // Allocate memory for g_lag
                state.g_lag.clear();
                    state.g_lag.push_back(X_Vector());
                    X::init(x,state.g_lag.back());

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

                // Set z to be ||x||/||e|| e.  In this way, ||z||=||x||.
                Z_Vector& zz=state.z.front();
                Z::id(zz);
                Z::scal(sqrt(X::innr(x,x))/sqrt(Z::innr(zz,zz)),zz);
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

            // Checks whether we have a valid real label.
            struct is_real : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( typename Unconstrained <Real,XX>::Restart
                            ::is_real()(name) ||
                        name == "mu" ||
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
                            ::is_x()(name) ||
                        name == "g_orig" ||
                        name == "g_schur" ||
                        name == "g_lag"
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
                xs.first.push_back("g_orig");
                xs.second.splice(xs.second.end(),state.g_orig);
                xs.first.push_back("g_schur");
                xs.second.splice(xs.second.end(),state.g_schur);
                xs.first.push_back("g_lag");
                xs.second.splice(xs.second.end(),state.g_lag);

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

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name0=="g_orig")
                        state.g_orig.splice(state.g_orig.end(),xs.second,x0);
                    else if(*name0=="g_schur")
                        state.g_schur.splice(state.g_schur.end(),xs.second,x0);
                    else if(*name0=="g_lag")
                        state.g_lag.splice(state.g_lag.end(),xs.second,x0);

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
            // Actual storage of the functions required
            struct t: public virtual Unconstrained <Real,XX>::Functions::t {
                // Inequality constraints 
                std::auto_ptr <VectorValuedFunction <Real,XX,ZZ> > h;
            };

            // A function with a modified Hessian used in inequality
            // constrained optimization
            struct InequalityModifiedFunction 
                : public peopt::ScalarValuedFunction <Real,XX>
            {
            private:
                //Underlying function.  
                std::auto_ptr <peopt::ScalarValuedFunction <Real,XX> > f;

                // Inequality constraint. For a reason why we have a reference
                // to an auto_ptr, see the discussion inside
                // Unconstrained::Functions::Merit
                const std::auto_ptr <peopt::VectorValuedFunction <Real,XX,ZZ> >&
                    h;

                // Inequality Lagrange multiplier
                const Z_Vector& z;

                // Interior point parameter
                const Real& mu;

                // Centering paramter 
                const Real& sigma;

                // Inequality constraint evaluated at x
                const Z_Vector& h_x;

                // Original version of the gradient
                mutable X_Vector& g_orig;

                // Schur complement version of the gradient
                mutable X_Vector& g_schur;

                // True gradient of the Lagrangian 
                mutable X_Vector& g_lag;

                // This forces derived classes to call the constructor that
                // depends on the state
                InequalityModifiedFunction() {}

            public:
                // The constructor determines whether we really need to build
                // a Hessian-vector product or if we use an internal
                // approximation
                InequalityModifiedFunction(
                    const Messaging& msg,
                    typename State::t& state,
                    typename Functions::t& fns
                ) : f(fns.f), h(fns.h), z(state.z.front()), mu(state.mu),
                    sigma(state.sigma), h_x(state.h_x.front()),
                    g_orig(state.g_orig.front()),
                    g_schur(state.g_schur.front()),
                    g_lag(state.g_lag.front())
                { }

                // <- f(x) 
                Real operator () (const X_Vector& x) const {
                    return (*f)(x);
                }

                // g = grad f(x) - sigma mu h'(x)* (inv(L(h(x))) e)
                void grad(const X_Vector& x,X_Vector& g) const {
                    // Create work elements for accumulating the
                    // interior point pieces
                    Z_Vector z_tmp1; Z::init(z,z_tmp1);
                    Z_Vector z_tmp2; Z::init(z,z_tmp2);
                    X_Vector x_tmp1; X::init(x,x_tmp1);

                    // g <- grad f(x)
                    f->grad(x,g);

                    // g_orig <- grad f(x)
                    X::copy(g,g_orig);

                    // Calculate the Schur complement version of the gradient

                    // g_schur <- grad f(x)
                    X::copy(g,g_schur);

                    // z_tmp1 <- e
                    Z::id(z_tmp1);

                    // z_tmp2 <- inv(L(h(x))) e 
                    Z::linv(h_x,z_tmp1,z_tmp2);

                    // x_tmp1 <- h'(x)* (inv(L(h(x))) e)
                    h->ps(x,z_tmp2,x_tmp1);

                    // g_schur <- grad f(x) - sigma mu h'(x)* (inv(L(h(x))) e)
                    X::axpy(-sigma*mu,x_tmp1,g_schur);

                    // Calculate the true gradient of the Lagrangian

                    // g_lag <- grad f(x)
                    X::copy(g,g_lag);

                    // x_tmp1 <- h'(x)*z
                    h->ps(x,z,x_tmp1);

                    // g_lag <- grad f(x) - h'(x)*z
                    X::axpy(-Real(1.0),x_tmp1,g_lag);

                    // Finally, return the Schur complement gradient
                    X::copy(g_schur,g);
                }

                // Adds on the interior point piece to the Hessian. 
                //
                // Primal-dual
                // H_dx
                // = hess f(x)dx + h'(x)* (inv(L(h(x))) (z o h'(x) dx)
                //                        + z o (inv L(h(x)) h'(x) dx))/2
                // ~= hess f(x) dx + h'(x)* inv(L(h(x))) (z o (h'(x) dx))
                //  = hess f(x) dx + h'(x)* inv(L(h(x))) L(z) h'(x) dx
                //
                // Log-barrier
                // Replace z with sigma mu inv(L(h(x)) e.  This comes from
                // the equation h(x) o z = sigma mu e.
                // 
                // Note, we symmetrize the interior point piece added since,
                // for SDP, this piece is likely not symmetric.
                virtual void hessvec(
                    const X_Vector& x,
                    const X_Vector& dx,
                    X_Vector& H_dx 
                ) const {
                    // Create work elements for accumulating the
                    // interior point pieces
                    Z_Vector z_tmp1; Z::init(z,z_tmp1);
                    Z_Vector z_tmp2; Z::init(z,z_tmp2);
                    Z_Vector linv_hx_z_prod_hpx; Z::init(z,linv_hx_z_prod_hpx);
                    X_Vector hess_mod; X::init(x,hess_mod);

                    // H_dx <- hess f(x) dx
                    f->hessvec(x,dx,H_dx);

                    // z_tmp1 <- h'(x) dx
                    h->p(x,dx,z_tmp1);

                    // z_tmp2 <- z o h'(x) dx
                    Z::prod(z,z_tmp1,z_tmp2);

                    // linv_hx_z_prod_hpx <- inv(L(h(x))) (z o h'(x) dx) 
                    Z::linv(h_x,z_tmp2,linv_hx_z_prod_hpx);
                    
                    // z_tmp2 <- inv L(h(x)) h'(x) dx
                    Z::linv(h_x,z_tmp1,z_tmp2);
                    
                    // z_tmp1 <- z o (inv L(h(x)) h'(x) dx)
                    Z::prod(z,z_tmp2,z_tmp1);

                    // z_tmp1 <- (inv(L(h(x))) (z o h'(x) dx)
                    //               + z o (inv L(h(x)) h'(x) dx))/2
                    Z::axpy(Real(1.),linv_hx_z_prod_hpx,z_tmp1);
                    Z::scal(Real(.5),z_tmp1);

                    // hess_mod <- h'(x)* (inv(L(h(x))) (z o h'(x) dx)
                    //                    + z o (inv L(h(x)) h'(x) dx))/2
                    h->ps(x,z_tmp1,hess_mod);

                    // H_dx 
                    //  = hess f(x) dx + h'(x)* (inv(L(h(x))) (z o h'(x) dx)
                    //                          + z o (inv L(h(x)) h'(x) dx))/2
                    X::axpy(Real(1.),hess_mod,H_dx);
                }
            };

            // A log-barrier merit function 
            struct Merit : public peopt::ScalarValuedFunction <Real,XX> {
            private:
                // Messaging tool
                const Messaging& msg;
                
                // Underlying function.  This takes control of the memory
                std::auto_ptr <peopt::ScalarValuedFunction <Real,XX> >
                    f_merit;

                // Inequality constraint.
                const std::auto_ptr <peopt::VectorValuedFunction <Real,XX,ZZ> >&
                    h;

                // Interior point parameter
                const Real& mu;

                // Interior point centering parameter
                const Real& sigma;

                // Inequality Lagrange multiplier
                const Z_Vector& z;
                
                // Make sure we call the constructor below
                Merit () {}
            public:
                Merit(
                    const Messaging& msg_,
                    const typename State::t& state,
                    typename Functions::t& fns
                ) : msg(msg_), f_merit(fns.f_merit), h(fns.h), mu(state.mu),
                    sigma(state.sigma), z(state.z.back()) {}

                // <- f_merit(x) - sigma mu barr(h(x))
                Real operator () (const X_Vector& x) const {
                    // Calculate h(x).  Funny enough, we don't use the cached
                    // value of h(x) here since we evaluate this function at
                    // points such as x+s.  Be careful of this.
                    Z_Vector h_x; Z::init(z,h_x);
                    (*h)(x,h_x);

                    // Return f_merit(x) - mu barr(h(x))
                    return (*f_merit)(x) - sigma*mu * Z::barr(h_x); 
                }

                // Throw an error 
                void grad(const X_Vector& x,X_Vector& g) const {
                    msg.error("The gradient of the merit function is "
                        "undefined.");
                }

                // Throw an error 
                void hessvec(
                    const X_Vector& x,
                    const X_Vector& dx,
                    X_Vector& H_dx
                ) const {
                    msg.error("The Hessian of the merit function is "
                        "undefined.");
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
                // Modify the objective 
                fns.f.reset(new InequalityModifiedFunction(msg,state,fns));

                // Create the merit function
                fns.f_merit.reset(new Merit(msg,state,fns));

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

            // Gets the header for the state information
            static void getStateHeader_(
                const typename State::t& state,
                std::list <std::string>& out
            ) {

                // Interior point pieces of the state 
                out.push_back(atos <> ("mu"));
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
                const typename State::t& state,
                const bool blank,
                std::list <std::string>& out
            ) {

                // Create some shortcuts
                const Real& mu=state.mu; 

                // Get a iterator to the last element prior to inserting
                // elements
                std::list <std::string>::iterator prior=out.end(); prior--;

                // Interior point information
                // Estimate the interior-point parameter
                out.push_back(atos <> (mu));

                // If we needed to do blank insertions, overwrite the elements
                // with spaces 
                if(blank)
                    for(std::list <std::string>::iterator x=++prior;
                        x!=out.end();
                        x++
                    )
                        (*x)="          ";
            }

            // Combines all of the state information
            static void getState(
                const typename State::t& state,
                const bool blank,
                const bool noiter,
                std::list <std::string>& out
            ) {
                Unconstrained <Real,XX>::Printer
                    ::getState_(state,blank,noiter,out);
                InequalityConstrained <Real,XX,ZZ>::Printer
                    ::getState_(state,blank,out);
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
            // z = inv L(h(x)) (-z o h'(x)s + sigma mu e)
            static void findInequalityMultiplierLinked(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& h_x=state.h_x.front();
                const X_Vector& x=state.x.front();
                const X_Vector& s=state.s.front();
                const Real& mu=state.mu;
                const Real& sigma=state.sigma;
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                Z_Vector& z=state.z.front();

                // z_tmp1 <- h'(x)s
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h.p(x,s,z_tmp1);

                // z_tmp2 <- z o h'(x)s
                Z_Vector z_tmp2; Z::init(z,z_tmp2);
                Z::prod(z,z_tmp1,z_tmp2);

                // z_tmp2 <- - z o h'(x)s
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 <- e
                Z::id(z_tmp1);

                // z_tmp2 <- -z o h'(x)s + sigma mu e
                Z::axpy(sigma*mu,z_tmp1,z_tmp2);

                // z <- inv L(h(x)) (-z o h'(x)s + sigma mu e)
                Z::linv(h_x,z_tmp2,z);
            }

            // Finds the new inequality Lagrange multiplier
            // z = sigma mu inv L(h(x)) e 
            static void findInequalityMultiplierLogBarrier(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& h_x=state.h_x.front();
                const Real& mu=state.mu;
                const Real& sigma=state.sigma;
                Z_Vector& z=state.z.front();

                // z_tmp1 <- e 
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                Z::id(z_tmp1);

                // z <- inv(L(h(x))) e 
                Z::linv(h_x,z_tmp1,z);

                // z <- sigma mu inv(L(h(x))) e 
                Z::scal(sigma*mu,z);
            }


            // Finds the new inequality Lagrange multiplier step
            // dz = -z + inv L(h(x)) (-z o h'(x)s + sigma mu e)
            static void findInequalityMultiplierStep(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();
                const X_Vector& x=state.x.front();
                const X_Vector& s=state.s.front();
                const Real& mu=state.mu;
                const Real& sigma=state.sigma;
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                Z_Vector& dz=state.dz.front();

                // z_tmp1 <- h'(x)s
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h.p(x,s,z_tmp1);

                // z_tmp2 <- z o h'(x)s
                Z_Vector z_tmp2; Z::init(z,z_tmp2);
                Z::prod(z,z_tmp1,z_tmp2);

                // z_tmp2 <- - z o h'(x)s
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 <- e
                Z::id(z_tmp1);

                // z_tmp2 <- -z o h'(x)s + sigma mu e
                Z::axpy(sigma*mu,z_tmp1,z_tmp2);

                // dz <- inv L(h(x)) (-z o h'(x)s + sigma mu e)
                Z::linv(h_x,z_tmp2,dz);

                // dz <- -z + inv L(h(x)) (-z o h'(x)s + sigma mu e)
                Z::axpy(Real(-1.),z,dz);
            }

            // Estimates the interior point parameter with the formula
            // mu = <z,h(x)>/m
            static Real estimateInteriorPointParameter(
                const typename Functions::t& fns,
                const typename State::t& state
            ) {
                // Create some shortcuts
                const Z_Vector& z=state.z.front();
                const Z_Vector& h_x=state.h_x.front();

                // Determine the scaling factor for the interior-
                // point parameter estimate
                Z_Vector z_tmp; Z::init(z,z_tmp);
                Z::id(z_tmp);
                Real m = Z::innr(z_tmp,z_tmp);

                // Estimate the interior-point parameter
                return Z::innr(z,h_x) / m;
            }
           
            // Adjust the stopping conditions unless mu < mu_typ*eps_mu 
            static void adjustStoppingConditions(
                const typename Functions::t& fns,
                typename State::t& state
            ) {
                // Create some shortcuts
                const Real& mu=state.mu;
                const Real& mu_typ=state.mu_typ;
                const Real& eps_mu=state.eps_mu;
                StoppingCondition::t& opt_stop=state.opt_stop;
                
                // Prevent convergence unless mu has been reduced to
                // eps_mu times mu_typ.
                if( opt_stop==StoppingCondition::RelativeGradientSmall &&
                    !(mu < mu_typ*eps_mu) 
                )
                    opt_stop=StoppingCondition::NotConverged;

                // If the interior point paramter is negative, exit
                if(mu < Real(0.))
                    opt_stop=StoppingCondition::InteriorPointInstability;
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
                X_Vector& s=state.s.front();
                Z_Vector& dz=state.dz.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector ss; X::init(x,ss);
                X::copy(s,ss);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,ss);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+s
                X_Vector x_tmp1; X::init(x,x_tmp1);
                X::copy(x,x_tmp1);
                X::axpy(Real(1.),ss,x_tmp1);

                // z_tmp1=h(x+s)
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h(x_tmp1,z_tmp1);

                // z_tmp1=h(x+s)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+s)-h(x)) + h(x) >=0
                Real alpha_x=Z::srch(z_tmp1,h_x);

                // Determine how far we can go in the dual variable 

                // Find the largest alpha such that
                // alpha dz + z >=0
                Real alpha_z=Z::srch(dz,z); 

                // Figure out how much to shorten the steps, if at all
                Real beta_x = alpha_x<Real(0.) || alpha_x*gamma>Real(1.)
                    ? Real(1.) : alpha_x*gamma;
                Real beta_z = alpha_z<Real(0.) || alpha_z*gamma>Real(1.)
                    ? Real(1.) : alpha_z*gamma;

                // Shorten the inequality multiplier step
                Z::scal(beta_z,dz);

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(beta_x,s);

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
                const Real& sigma=state.sigma;
                const AlgorithmClass::t& algorithm_class =state.algorithm_class;
                const Z_Vector& z=state.z.front();
                const X_Vector& x=state.x.front();
                const Z_Vector& h_x=state.h_x.front();
                const VectorValuedFunction <Real,XX,ZZ>& h=*(fns.h);
                X_Vector& s=state.s.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector ss; X::init(x,ss);
                X::copy(s,ss);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,ss);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+s
                X_Vector x_tmp1; X::init(x,x_tmp1);
                X::copy(x,x_tmp1);
                X::axpy(Real(1.),ss,x_tmp1);

                // z_tmp1=h(x+s)
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h(x_tmp1,z_tmp1);

                // z_tmp2=h(x+s)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+s)-h(x)) + h(x) >=0
                Real alpha1=Z::srch(z_tmp1,h_x);

                // Determine how far we can go in the dual variable

                // z_tmp1=h'(x)s
                h.p(x,ss,z_tmp1);
                
                // z_tmp2 = z o h'(x)s
                Z_Vector z_tmp2; Z::init(z,z_tmp2);
                Z::prod(z,z_tmp1,z_tmp2);

                // z_tmp2 = - z o h'(x)s
                Z::scal(Real(-1.),z_tmp2);

                // z_tmp1 = e
                Z::id(z_tmp1);

                // z_tmp1 = sigma mu e
                Z::scal(sigma*mu,z_tmp1);

                // Find the largest alpha such that
                // alpha (- z o h'(x)s) + sigma mu e >=0
                Real alpha2=Z::srch(z_tmp2,z_tmp1);

                // Determine the farthest we can go in both variables
                Real alpha0;

                // Only the dual step is restrictive
                if(alpha1 < Real(0.) && alpha2 >= Real(0.))
                    alpha0 = alpha2;

                // Only the primal step is restrictive
                else if(alpha1 >= Real(0.) && alpha2 < Real(0.))
                    alpha0 = alpha1;

                // Neither step is restrictive
                else if(alpha1 < Real(0.) && alpha2 < Real(0.))
                    alpha0 = Real(-1.);

                // Both steps are restrictive
                else
                    alpha0 = alpha1 < alpha2 ? alpha1 : alpha2;
                    
                // Next, determine if we need to back off from the
                // boundary or leave the step unchanged.
                alpha0 =
                    alpha0 == Real(-1.) || alpha0*gamma > Real(1.) ?
                        // If we're unrestricted or backing off still puts
                        // us farther than we want to go
                        Real(1.) :
                        // If we're restricted, back off from the boundary
                        alpha0*gamma;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(alpha0,s);

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
                X_Vector& s=state.s.front();
                Real& alpha=state.alpha;

                // Create a fake step.  In the case of a trust-region
                // method this is just the step.  In the case of
                // a line-search method this is 2 alpha s.  This represents
                // the farthest either method will attempt to step.
                X_Vector ss; X::init(x,ss);
                X::copy(s,ss);
                if(algorithm_class==AlgorithmClass::LineSearch)
                    X::scal(Real(2.)*alpha,ss);
                
                // Determine how far we can go in the primal variable
                
                // x_tmp1=x+s
                X_Vector x_tmp1; X::init(x,x_tmp1);
                X::copy(x,x_tmp1);
                X::axpy(Real(1.),ss,x_tmp1);

                // z_tmp1=h(x+s)
                Z_Vector z_tmp1; Z::init(z,z_tmp1);
                h(x_tmp1,z_tmp1);

                // z_tmp1=h(x+s)-h(x)
                Z::axpy(Real(-1.),h_x,z_tmp1);

                // Find the largest alpha such that
                // alpha (h(x+s)-h(x)) + h(x) >=0
                Real alpha_x=Z::srch(z_tmp1,h_x);

                // Figure out how much to shorten the steps, if at all
                Real beta_x = alpha_x<Real(0.) || alpha_x*gamma>Real(1.)
                    ? Real(1.) : alpha_x*gamma;

                // If we're doing a trust-region method, shorten the
                // step length accordingly
                if(algorithm_class==AlgorithmClass::TrustRegion) 

                    // Shorten the step
                    X::scal(beta_x,s);

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
                    const X_Vector& g_orig=state.g_orig.front();
                    const X_Vector& g_lag=state.g_lag.front();
                    const X_Vector& g_schur=state.g_schur.front();
                    const InteriorPointMethod::t& ipm=state.ipm;
                    const CentralityStrategy::t& cstrat=state.cstrat;
                    const unsigned int& iter=state.iter;
                    X_Vector& g=state.g.front();
                    Z_Vector& z=state.z.front();
                    Z_Vector& h_x=state.h_x.front();
                    Z_Vector& dz=state.dz.front();
                    Real& mu = state.mu;
                    Real& sigma = state.sigma;
                    Real& merit_xps= state.merit_xps;
                    Real& norm_g = state.norm_g;
                    const ScalarValuedFunction<Real,XX>& f_merit=*(fns.f_merit);

                    // Call the user define manipulator
                    smanip(fns,state,loc);

                    switch(loc){
                    case OptimizationLocation::BeforeInitialFuncAndGrad:
                        // Do predictor corrector if need be
                        if(cstrat==CentralityStrategy::PredictorCorrector)
                            sigma = iter % 2 ? Real(0.) : Real(1.);

                        // In a log-barrier method, find the initial Lagrange
                        // multiplier.
                        if(ipm==InteriorPointMethod::LogBarrier)
                            findInequalityMultiplierLogBarrier(fns,state);
                        break;

                    // After we reject a step, get the gradient of the
                    // Lagrangian for the output.  In addition, make sure
                    // we specify that we take a zero step in the Lagrange
                    // multiplier.  This is important in case we exit early
                    // due to small steps.
                    case OptimizationLocation::AfterRejectedTrustRegion:
                    case OptimizationLocation::AfterRejectedLineSearch:
                        X::copy(g_lag,g);
                        norm_g=sqrt(X::innr(g_lag,g_lag));

                        Z::zero(dz);
                        break;

                    // Make sure we use the Schur complement gradient when
                    // we get the step
                    case OptimizationLocation::BeforeGetStep:
                        X::copy(g_schur,g);
                        norm_g=sqrt(X::innr(g_schur,g_schur));
                        break;
                    
                    case OptimizationLocation::BeforeStep:
                        // Do predictor corrector if need be
                        if(cstrat==CentralityStrategy::PredictorCorrector)
                            // We have a iter+1 here since we haven't updated
                            // our iteration yet and we need to set sigma
                            // for the next iteration.
                            sigma = (iter+1) % 2 ? Real(0.) : Real(1.);

                        // Find the new inequality multiplier or step
                        switch(ipm){
                        case InteriorPointMethod::PrimalDual:
                            Z::axpy(Real(1.),dz,z);
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
                        mu = estimateInteriorPointParameter(fns,state);

                        // Update the inequality multiplier in a log-barrier
                        // method
                        if(ipm==InteriorPointMethod::LogBarrier)
                            findInequalityMultiplierLogBarrier(fns,state);

                        // Calculate the new merit function based on this
                        // interior point paramter
                        merit_xps=f_merit(x);
                        break; 

                    // Copy in the original gradient for the quasi-Newton method
                    // or a CG-style algorithm.
                    case OptimizationLocation::BeforeSaveOld:
                    case OptimizationLocation::BeforeQuasi:
                        X::copy(g_orig,g);
                        norm_g=sqrt(X::innr(g_orig,g_orig));
                        break;

                    // Copy in the gradient of the Lagrangian for the
                    // stopping conditions.
                    case OptimizationLocation::AfterQuasi:
                        X::copy(g_lag,g);
                        norm_g=sqrt(X::innr(g_lag,g_lag));
                        break;

                    // Adjust the interior point parameter and insure that
                    // we do not converge unless the interior point parameter
                    // is small.
                    case OptimizationLocation::EndOfOptimizationIteration:
                        adjustStoppingConditions(fns,state);
                        break;

                    // Adjust our step or potential step to preserve positivity
                    case OptimizationLocation::BeforeLineSearch:
                    case OptimizationLocation::BeforeActualVersusPredicted:
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

                // Add the interior point pieces to the state manipulator
                InteriorPointManipulator <InequalityConstrained <Real,XX,ZZ> >
                    ipmanip(smanip);
                
                // Adds the output pieces to the state manipulator 
                DiagnosticManipulator <InequalityConstrained <Real,XX,ZZ> >
                    dmanip(ipmanip,msg);

                // Insures that we can interact with unconstrained code
                ConversionManipulator
                    <InequalityConstrained<Real,XX,ZZ>,Unconstrained <Real,XX> >
                    cmanip(dmanip);
                
                // Initialize any remaining functions required for optimization 
                Functions::init(msg,state,fns);

                // Initialize the value h(x)
                (*fns.h)(state.x.back(),state.h_x.back());

                // Estimate the interior point parameter
                state.mu=estimateInteriorPointParameter(fns,state);

                // Set the typical value for mu
                state.mu_typ=state.mu;

                // Set the trust-region scaling
                fns.TR_op.reset(new TrustRegionScaling(fns,state));
                
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
                            std::list <unsigned int> > Nats;
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

            // The actual internal state of the optimization
            struct t: 
                public EqualityConstrained <Real,XX,YY>::State::t,
                public InequalityConstrained <Real,XX,ZZ>::State::t
            {
                // Initialization constructors
                t() {
                    Constrained <Real,XX,YY,ZZ>::State::init_params(*this);
                }
                t(const X_Vector& x,const Y_Vector& y,const Z_Vector& z) {
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
                EqualityConstrained <Real,XX,YY>::State::init_vectors_(state,y);
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
                    ::Restart::stateToVectors(state,ys);
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
                    ::Restart::vectorsToState(state,ys);
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

            // Actual storage of the functions required
            struct t: 
                public EqualityConstrained <Real,XX,YY>::Functions::t,
                public InequalityConstrained <Real,XX,ZZ>::Functions::t
            {};

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

    };
}
#endif
