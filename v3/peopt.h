#ifndef PEOPT_H
#define PEOPT_H

#include<list>
#include<vector>
#include<map>
#include<limits>
#include<cmath>
#include<sstream>
#include<iomanip>
#include<memory>

namespace peopt{

    // A trick to determine if two types are the same
    template <typename T1,typename T2>
    struct is_same{
    private:
        static const bool same=false;
    public:
        static bool eval() {return same;}
    };
    template <typename T>
    struct is_same <T,T>{
    private:
        static const bool same=true;
    public:
        static bool eval() {return same;}
    };

    // Defines a Hilbert space 
    template <typename Vector_,typename Real_>
    struct HilbertSpace {
        // Define the vector space and the field
        typedef Vector_ Vector;
        typedef Real_ Real;

        // Create an empty, uninitialized vector
        virtual Vector create() = 0;

        // Memory allocation and size setting
        virtual void init(const Vector& x,Vector& y) = 0;

        // y <- x (Shallow.  No memory allocation.)
        virtual void copy(const Vector& x,Vector& y) = 0;

        // x <- alpha * x
        virtual void scal(const Real& alpha,Vector& x) = 0;

        // x <- 0.  Part of the reason we have this function and not just
        // use scal is that if the elements of x become NaN, then the
        // scaling operation may be undefined.  This is a hard set, which
        // should always be safe.
        virtual void zero(Vector& x) = 0;

        // y <- alpha * x + y
        virtual void axpy(const Real& alpha,const Vector& x,Vector& y) = 0;

        // innr <- <x,y>
        virtual Real innr(const Vector& x,const Vector& y) = 0;

        // norm <- ||x||
        Real norm(const Vector& x) {
            return sqrt(innr(x,x));
        }

        // Deallocates memory required by the vector space
        virtual ~HilbertSpace() {}
    };
    
    // Augments a Hilbert space with a Euclidean-Jordean algebra that defines
    // a cone
    template <typename Vector_,typename Real_>
    struct Cone : public HilbertSpace <Vector_,Real_> {
        // Define the vector space and the field
        typedef Vector_ Vector;
        typedef Real_ Real;

        // Jordan product, z <- x o y
        virtual void prod(
            const Vector& x,
            const Vector& y,
            Vector& z
        ) const = 0;

        // Identity element, x <- e such that x o e = x
        virtual void id(Vector& x) const = 0;

        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
        virtual void linv(
            const Vector& x,
            const Vector& y,
            Vector& z
        ) const = 0;
    };

    // Defines how we output messages to the user
    struct Messaging {
        // Defines the current print level
        unsigned int plevel;

        // Sets the default messaging level to 1
        Messaging() : plevel(1) {};

        // Changing the default messaging level during construction
        Messaging(unsigned int plevel_) : plevel(plevel_) {};

        // Prints a message
        virtual void print(const std::string msg,const unsigned int level)
        const {
            if(level <= plevel) std::cout << msg << std::endl;
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
            LineSearch              // Line-search algorithms
        };

        // Converts the algorithm class to a string
        std::string to_string(t algorithm_class){
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
        t from_string(std::string algorithm_class){
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
    }

    // Reasons why we stop the algorithm
    namespace StoppingCondition{
        enum t{
            NotConverged,            // Algorithm did not converge
            RelativeGradientSmall,   // Relative gradient was sufficiently small
            RelativeStepSmall,       // Relative change in the step is small
            MaxItersExceeded,        // Maximum number of iterations exceeded
            External                 // Some external stopping condition 
        };

        // Converts the stopping condition to a string 
        std::string to_string(t opt_stop){
            switch(opt_stop){
            case NotConverged:
                return "NotConverged";
            case RelativeGradientSmall:
                return "RelativeGradientSmall";
            case RelativeStepSmall:
                return "RelativeStepSmall";
            case MaxItersExceeded:
                return "MaxItersExceeded";
            case External:
                return "External";
            default:
                throw;
            }
        }

        // Converts a string to a stopping condition
        t from_string(std::string opt_stop){
            if(opt_stop=="NotConverged")
                return NotConverged;
            else if(opt_stop=="RelativeGradientSmall")
                return RelativeGradientSmall;
            else if(opt_stop=="RelativeStepSmall")
                return RelativeStepSmall;
            else if(opt_stop=="MaxItersExceeded")
                return MaxItersExceeded;
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
                    name=="External"
                )
                    return true;
                else
                    return false;
            }
        };
    }

    // Reasons we stop the Krylov method
    namespace KrylovStop{
        enum t{
            NegativeCurvature,        // Negative curvature detected
            RelativeErrorSmall,       // Relative error is small
            MaxItersExceeded,         // Maximum number of iterations exceeded
            TrustRegionViolated       // Trust-region radius violated
        };

        // Converts the Krylov stopping condition to a string 
        std::string to_string(t krylov_stop){
            switch(krylov_stop){
            case NegativeCurvature:
                return "NegativeCurvature";
            case RelativeErrorSmall:
                return "RelativeErrorSmall";
            case MaxItersExceeded:
                return "MaxItersExceeded";
            case TrustRegionViolated:
                return "TrustRegionViolated";
            default:
                throw;
            }
        }
        
        // Converts a string to a Krylov stopping condition
        t from_string(std::string krylov_stop){
            if(krylov_stop=="NegativeCurvature")
                return NegativeCurvature;
            else if(krylov_stop=="RelativeErrorSmall")
                return RelativeErrorSmall;
            else if(krylov_stop=="MaxItersExceeded")
                return MaxItersExceeded;
            else if(krylov_stop=="TrustRegionViolated")
                return TrustRegionViolated;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="NegativeCurvature" ||
                    name=="RelativeErrorSmall" ||
                    name=="MaxItersExceeded" ||
                    name=="TrustRegionViolated" 
                )
                    return true;
                else
                    return false;
            }
        };
    }

    // Various operators for both Hessian approximations and preconditioners
    namespace Operators{
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
        std::string to_string(t op){
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
        t from_string(std::string op){
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
    }

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
        std::string to_string(t dir){
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
        t from_string(std::string dir){
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
    }

    namespace LineSearchKind{
        enum t{
            Brents,           // Brent's minimization
            GoldenSection,    // Golden-section search 
            BackTracking,     // BackTracking search 
            TwoPointA,        // Barzilai and Borwein's method A
            TwoPointB         // Barzilai and Borwein's method B
        };
            
        // Converts the line-search kind to a string 
        std::string to_string(t kind){
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
        t from_string(std::string kind){
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
    }
    
    namespace OptimizationLocation{
        enum t{
            // Occurs after we take the optimization step u+s, but before
            // we calculate the gradient based on this new step.  In addition,
            // after this point we set the objective value, obj_u, to be
            // obj_ups.
            AfterStepBeforeGradient,

            // This occurs last in the optimization loop.  At this point,
            // we have already incremented our optimization iteration and
            // checked our stopping condition
            EndOfOptimizationIteration  
        };
            
        // Converts the line-search kind to a string 
        std::string to_string(t loc){
            switch(loc){
            case AfterStepBeforeGradient:
                return "AfterStepBeforeGradient";
            case EndOfOptimizationIteration:
                return "EndOfOptimizationIteration";
            default:
                throw;
            }
        }
        
        // Converts a string to a line-search kind 
        t from_string(std::string loc){
            if(loc=="AfterStepBeforeGradient")
                return AfterStepBeforeGradient; 
            else if(loc=="EndOfOptimizationIteration")
                return EndOfOptimizationIteration; 
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="AfterStepBeforeGradient" ||
                    name=="EndOfOptimizationIteration"
                )
                    return true;
                else
                    return false;
            }
        };
    }

    // Different nonlinear-CG directions 
    namespace NonlinearCGDirections {
        enum t{
            HestenesStiefel,        
            PolakRibiere,           
            FletcherReeves          
        };
    }

    // The structures in this namespace represent the internal state of the
    // optimization algorithm.
    namespace State {

        // State of an unconstrained optimization problem of the form
        // 
        // min_{x \in X} f(x)
        //
        // where f : X -> R
        template <typename Real,template <typename> class XX> 
        struct Unconstrained {
        public:
            // Create some shortcuts for some type names
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;

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

            // Current objective value
            Real obj_x;

            // Objective value at the trial step
            Real obj_xps;
            
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
            
            // ------------- CONSTRUCTORS ------------ 

            // Initialize the state without setting up any variables. 
            Unconstrained() {
                init_params();
            };

            // Initialize the state for unconstrained optimization.  
            Unconstrained(const X_Vector& x) {
                init_params();
                init_vectors(x);
            }

        protected:
            // This initializes all the variables required for unconstrained
            // optimization.  These variables are also require for constrained
            // optimization
            void init_vectors(const X_Vector& x_) {
                x.push_back(X::create()); X::init(x_,x.back());
                    X::copy(x_,x.back());
                g.push_back(X::create()); X::init(x_,g.back());
                s.push_back(X::create()); X::init(x_,s.back()); 
                x_old.push_back(X::create()); X::init(x_,x_old.back()); 
                g_old.push_back(X::create()); X::init(x_,g_old.back()); 
                s_old.push_back(X::create()); X::init(x_,s_old.back()); 
            }

            // This sets all of the parameters possible that don't require
            // special memory allocation such as variables.
            void init_params(){
                eps_g=Real(1e-6);
                eps_s=Real(1e-6);
                stored_history=0;
                history_reset=5;
                iter=1;
                iter_max=10;
                opt_stop=StoppingCondition::NotConverged;
                krylov_iter=1;
                krylov_iter_max=10;
                krylov_iter_total=0;
                krylov_stop=KrylovStop::RelativeErrorSmall;
                krylov_rel_err=Real(std::numeric_limits<double>::quiet_NaN());
                eps_krylov=Real(1e-2);
                algorithm_class=AlgorithmClass::TrustRegion;
                Minv_type=Operators::Identity;
                H_type=Operators::Identity;
                norm_g=Real(std::numeric_limits<double>::quiet_NaN());
                norm_gtyp=Real(std::numeric_limits<double>::quiet_NaN());
                norm_s=Real(std::numeric_limits<double>::quiet_NaN());
                norm_styp=Real(std::numeric_limits<double>::quiet_NaN());
                obj_x=Real(std::numeric_limits<double>::quiet_NaN());
                obj_xps=Real(std::numeric_limits<double>::quiet_NaN());
                delta=Real(100.);
                delta_max=Real(100.);
                eta1=Real(.1);
                eta2=Real(.9);
                rho=Real(0.);
                rejected_trustregion=0;
                alpha=1.;
                linesearch_iter=0;
                linesearch_iter_max=5;
                linesearch_iter_total=0;
                eps_ls=Real(1e-2);
                dir=LineSearchDirection::SteepestDescent;
                kind=LineSearchKind::GoldenSection;
            }

        public:
            
            // --------- VERIFYING PARAMETERS -------- 

            // Check that we have a valid set of parameters
            void check(const Messaging& msg){

                // Check that the tolerance for the gradient stopping condition
                // is positive
                if(eps_g <= 0) {
                    std::stringstream ss;
                    ss << "The tolerance for the gradient stopping condition "
                        "must be positive: eps_g = " << eps_g;
                    msg.error(ss.str());
                }
            
                // Check that the tolerance for the step length stopping
                // condition is positive
                if(eps_s <= 0) {
                    std::stringstream ss;
                    ss << "The tolerance for the step length stopping "
                        "condition must be positive: eps_s = " << eps_s; 
                    msg.error(ss.str());
                }
        
                // Check that the current iteration is positive
                if(iter <= 0) {
                    std::stringstream ss;
                    ss << "The current optimization iteration must be "
                        "positive: iter = " << iter;
                    msg.error(ss.str());
                }

                // Check that the maximum iteration is positive
                if(iter_max <= 0) {
                    std::stringstream ss;
                    ss << "The maximum optimization iteration must be "
                        "positive: iter_max = " << iter_max;
                    msg.error(ss.str());
                }

                // Check that the current Krylov iteration is positive
                if(krylov_iter <= 0) {
                    std::stringstream ss;
                    ss << "The current Krlov iteration must be "
                        "positive: krylov_iter = " << krylov_iter;
                    msg.error(ss.str());
                }

                // Check that the maximum Krylov iteration is positive
                if(krylov_iter_max <= 0) {
                    std::stringstream ss;
                    ss << "The maximum Krylov iteration must be "
                        "positive: krylov_iter_max = " << krylov_iter_max;
                    msg.error(ss.str());
                }

                // Check that relative error in the Krylov method is nonnegative
                if(krylov_rel_err < 0) {
                    std::stringstream ss;
                    ss << "The relative error in the Krylov method must be "
                        "nonnegative: krylov_rel_err = " << krylov_rel_err;
                    msg.error(ss.str());
                }
                
                // Check that the stopping tolerance for the Krylov method is
                // positive
                if(eps_krylov <= 0) {
                    std::stringstream ss;
                    ss << "The tolerance for the Krylov method stopping "
                        "condition must be positive: eps_krylov = "<<eps_krylov;
                    msg.error(ss.str());
                }

                // Check that the norm of the gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                if(norm_g < 0 || (iter!=1 && norm_g!=norm_g)) {
                    std::stringstream ss;
                    ss << "The norm of the gradient must be nonnegative: "
                        "norm_g = " << norm_g; 
                    msg.error(ss.str());
                }

                // Check that the norm of a typical gradient is nonnegative or
                // if we're on the first iteration, we allow a NaN
                if(norm_gtyp < 0 || (iter!=1 && norm_gtyp!=norm_gtyp)) {
                    std::stringstream ss;
                    ss << "The norm of a typical gradient must be nonnegative: "
                        "norm_gtyp = " << norm_gtyp; 
                    msg.error(ss.str());
                }

                // Check that the norm of the trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                if(norm_s < 0 || (iter!=1 && norm_s!=norm_s)) {
                    std::stringstream ss;
                    ss << "The norm of the trial step must be nonnegative: "
                        "norm_s = " << norm_s; 
                    msg.error(ss.str());
                }

                // Check that the norm of a typical trial step is nonnegative or
                // if we're on the first iteration, we allow a NaN
                if(norm_styp < 0 || (iter!=1 && norm_styp!=norm_styp)) {
                    std::stringstream ss;
                    ss << "The norm of a typical trial step must be "
                        "nonnegative: norm_styp = " << norm_styp; 
                    msg.error(ss.str());
                }

                // Check that the objective value isn't a NaN past iteration 1
                if(iter!=1 && obj_x!=obj_x) {
                    std::stringstream ss;
                    ss << "The objective value must be a number: obj_x = "
                        << obj_x;
                    msg.error(ss.str());
                }

                // Check that the objective at a trial step isn't a NaN past
                // iteration 1
                if(iter!=1 && obj_xps!=obj_xps) {
                    std::stringstream ss;
                    ss << "The objective value at the trial step must be a "
                        "number: obj_xps = " << obj_xps;
                    msg.error(ss.str());
                }

                // Check that the trust-region radius is positive
                if(delta<=0){
                    std::stringstream ss;
                    ss << "The trust-region radius must be positive: delta = "
                        << delta; 
                    msg.error(ss.str());
                }

                // Check that the maximum trust-region radius is positive
                if(delta_max<=0){
                    std::stringstream ss;
                    ss << "The maximum trust-region radius must be positive: "
                        "delta_max = " << delta_max; 
                    msg.error(ss.str());
                }

                // Check that the current trust-region radius is less than
                // or equal to the maximum trust-region radius
                if(delta > delta_max){
                    std::stringstream ss;
                    ss << "The trust-region radius must be less than or equal "
                        "to the maximum trust-region radius: delta = "
                        << delta << ", delta_max = " << delta_max;
                    msg.error(ss.str());
                }

                // Check that the predicted vs. actual reduction tolerance
                // is between 0 and 1
                if(eta1 < 0 || eta1 > 1){
                    std::stringstream ss;
                    ss << "The tolerance for whether or not we accept a "
                        "trust-region step must be between 0 and 1: eta1 = "
                        << eta1;
                    msg.error(ss.str());
                }
                
                // Check that the other predicted vs. actual reduction tolerance
                // is between 0 and 1
                if(eta2 < 0 || eta2 > 1){
                    std::stringstream ss;
                    ss << "The tolerance for whether or not we increase the "
                        "trust-region radius must be between 0 and 1: eta2 = "
                        << eta2;
                    msg.error(ss.str());
                }

                // Check that eta2 > eta1
                if(eta1 >= eta2) {
                    std::stringstream ss;
                    ss << "The trust-region tolerances for accepting steps "
                        "must satisfy the relationship that eta1 < eta2: "
                        "eta1 = " << eta1 << ", eta2 = " << eta2;
                    msg.error(ss.str());
                }

                // Check that the prediction versus actual reduction is
                // nonnegative 
                if(rho < 0) {
                    std::stringstream ss;
                    ss << "The predicted versus actual reduction must be "
                        "nonnegative: rho = " << rho;
                    msg.error(ss.str());
                }

                // Check that the line-search step length is positive 
                if(alpha <= 0) {
                    std::stringstream ss;
                    ss << "The line-search step length must be positive: "
                        "alpha = " << alpha;
                    msg.error(ss.str());
                }
                
                // Check that the stopping tolerance for the line-search
                // methods is positive
                if(eps_ls <= 0) {
                    std::stringstream ss;
                    ss << "The tolerance for the line-search stopping "
                        "condition must be positive: eps_ls = " << eps_ls;
                    msg.error(ss.str());
                }
            }

            // -------------- RESTARTING ------------- 

        public:
            // Create a series of types used for peering, capturing, and
            // releasing information in the state
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <unsigned int> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 

        protected:
            // ----- COPYING TO AND FROM STATE ------- 

            // Copy out all variables.
            void stateToVectors(X_Vectors& xs) {
                // Move the memory of all variables into the list 
                xs.first.push_back("x");
                xs.second.splice(xs.second.end(),x);
                xs.first.push_back("g");
                xs.second.splice(xs.second.end(),g);
                xs.first.push_back("s");
                xs.second.splice(xs.second.end(),s);
                xs.first.push_back("x_old");
                xs.second.splice(xs.second.end(),x_old);
                xs.first.push_back("g_old");
                xs.second.splice(xs.second.end(),g_old);
                xs.first.push_back("s_old");
                xs.second.splice(xs.second.end(),s_old);

                // Write out the quasi-Newton information with sequential names
                {int i=1;
                for(typename std::list <X_Vector>::iterator y=oldY.begin();
                    y!=oldY.end();
                    y=oldY.begin()
                ){
                    std::stringstream ss;
                    ss << "oldY_" << i;
                    xs.first.push_back(ss.str());
                    xs.second.splice(xs.second.end(),oldY,y);
                }}

                // Write out the quasi-Newton information with sequential names
                {int i=1;
                for(typename std::list <X_Vector>::iterator s=oldS.begin();
                    s!=oldS.end();
                    s=oldS.begin()
                ){
                    std::stringstream ss;
                    ss << "oldS_" << i;
                    xs.first.push_back(ss.str());
                    xs.second.splice(xs.second.end(),oldS,s);
                }}
            }
            
            // Copy out all non-variables.  This includes reals, naturals,
            // and parameters
            void stateToScalars(
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                
                // Copy in all the real numbers 
                reals.first.push_back("eps_g");
                reals.second.push_back(eps_g);
                reals.first.push_back("eps_s");
                reals.second.push_back(eps_s);
                reals.first.push_back("krylov_rel_err");
                reals.second.push_back(krylov_rel_err);
                reals.first.push_back("eps_krylov");
                reals.second.push_back(eps_krylov);
                reals.first.push_back("norm_g");
                reals.second.push_back(norm_g);
                reals.first.push_back("norm_gtyp");
                reals.second.push_back(norm_gtyp);
                reals.first.push_back("norm_s");
                reals.second.push_back(norm_s);
                reals.first.push_back("norm_styp");
                reals.second.push_back(norm_styp);
                reals.first.push_back("obj_x");
                reals.second.push_back(obj_x);
                reals.first.push_back("obj_xps");
                reals.second.push_back(obj_xps);
                reals.first.push_back("delta");
                reals.second.push_back(delta);
                reals.first.push_back("delta_max");
                reals.second.push_back(delta_max);
                reals.first.push_back("eta1");
                reals.second.push_back(eta1);
                reals.first.push_back("eta2");
                reals.second.push_back(eta2);
                reals.first.push_back("rho");
                reals.second.push_back(rho);
                reals.first.push_back("alpha");
                reals.second.push_back(alpha);
                reals.first.push_back("eps_ls");
                reals.second.push_back(eps_ls);

                // Copy in all the natural numbers
                nats.first.push_back("stored_history");
                nats.second.push_back(stored_history);
                nats.first.push_back("history_reset");
                nats.second.push_back(history_reset);
                nats.first.push_back("iter");
                nats.second.push_back(iter);
                nats.first.push_back("iter_max");
                nats.second.push_back(iter_max);
                nats.first.push_back("krylov_iter");
                nats.second.push_back(krylov_iter);
                nats.first.push_back("krylov_iter_max");
                nats.second.push_back(krylov_iter_max);
                nats.first.push_back("krylov_iter_total");
                nats.second.push_back(krylov_iter_total);
                nats.first.push_back("rejected_trustregion");
                nats.second.push_back(rejected_trustregion);
                nats.first.push_back("linesearch_iter");
                nats.second.push_back(linesearch_iter);
                nats.first.push_back("linesearch_iter_max");
                nats.second.push_back(linesearch_iter_max);
                nats.first.push_back("linesearch_iter_total");
                nats.second.push_back(linesearch_iter_total);

                // Copy in all the parameters
                params.first.push_back("algorithm_class");
                params.second.push_back(
                    AlgorithmClass::to_string(algorithm_class));
                params.first.push_back("opt_stop");
                params.second.push_back(StoppingCondition::to_string(opt_stop));
                params.first.push_back("krylov_stop");
                params.second.push_back(KrylovStop::to_string(krylov_stop));
                params.first.push_back("H_type");
                params.second.push_back(Operators::to_string(H_type));
                params.first.push_back("Minv_type");
                params.second.push_back(Operators::to_string(Minv_type));
                params.first.push_back("dir");
                params.second.push_back(LineSearchDirection::to_string(dir));
                params.first.push_back("kind");
                params.second.push_back(LineSearchKind::to_string(kind));
            }

            // Copy in all variables.  This assumes that the quasi-Newton
            // information is being read in order.
            void vectorsToState(X_Vectors& xs) {

                for(std::list <std::string>::iterator name=xs.first.begin();
                    name!=xs.first.end();
                    name++
                ) {
                    // Since we're using a splice operation, we slowly empty
                    // the variable list.  Hence, we always take the first
                    // element.
                    typename std::list <X_Vector>::iterator xx
                        =xs.second.begin();

                    // Determine which xxiable we're reading in and then splice
                    // it in the correct location
                    if(*name=="x") x.splice(x.end(),xs.second,xx);
                    else if(*name=="g") g.splice(g.end(),xs.second,xx);
                    else if(*name=="s") s.splice(s.end(),xs.second,xx); 
                    else if(*name=="x_old")
                        x_old.splice(x_old.end(),xs.second,xx);
                    else if(*name=="g_old")
                        g_old.splice(g_old.end(),xs.second,xx);
                    else if(*name=="s_old")
                        s_old.splice(s_old.end(),xs.second,xx);
                    else if(name->substr(0,5)=="oldY_")
                        oldY.splice(oldY.end(),xs.second,xx);
                    else if(name->substr(0,5)=="oldS_")
                        oldS.splice(oldS.end(),xs.second,xx);
                }
            }

            // Copy in all non-variables.  This includes reals, naturals,
            // and parameters
            void scalarsToState(
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
                    if(*name=="eps_g") eps_g=*real;
                    else if(*name=="eps_s") eps_s=*real;
                    else if(*name=="krylov_rel_err") krylov_rel_err=*real;
                    else if(*name=="eps_krylov") eps_krylov=*real;
                    else if(*name=="norm_g") norm_g=*real;
                    else if(*name=="norm_gtyp") norm_gtyp=*real;
                    else if(*name=="norm_s") norm_g=*real;
                    else if(*name=="norm_styp") norm_gtyp=*real;
                    else if(*name=="obj_x") obj_x=*real;
                    else if(*name=="obj_xps") obj_xps=*real;
                    else if(*name=="delta") delta=*real;
                    else if(*name=="delta_max") delta_max=*real;
                    else if(*name=="eta1") eta1=*real;
                    else if(*name=="eta2") eta2=*real;
                    else if(*name=="rho") rho=*real;
                    else if(*name=="alpha") alpha=*real;
                    else if(*name=="eps_ls") eps_ls=*real;
                }
            
                // Next, copy in any naturals
                std::list <unsigned int>::iterator nat=nats.second.begin();
                for(std::list <std::string>::iterator name=nats.first.begin();
                    name!=nats.first.end();
                    name++,nat++
                ){
                    if(*name=="stored_history") stored_history=*nat;
                    else if(*name=="history_reset") history_reset=*nat;
                    else if(*name=="iter") iter=*nat;
                    else if(*name=="iter_max") iter_max=*nat;
                    else if(*name=="krylov_iter") krylov_iter=*nat;
                    else if(*name=="krylov_iter_max") krylov_iter_max=*nat;
                    else if(*name=="krylov_iter_total") krylov_iter_total=*nat;
                    else if(*name=="rejected_trustregion")
                        rejected_trustregion=*nat;
                    else if(*name=="linesearch_iter") linesearch_iter=*nat;
                    else if(*name=="linesearch_iter_max")
                        linesearch_iter_max=*nat;
                    else if(*name=="linesearch_iter_total")
                        linesearch_iter_total=*nat;
                }
                    
                // Next, copy in any parameters 
                std::list <std::string>::iterator param=params.second.begin();
                for(std::list <std::string>::iterator name=params.first.begin();
                    name!=params.first.end();
                    name++,param++
                ){
                    if(*name=="algorithm_class")
                        algorithm_class=AlgorithmClass::from_string(*param);
                    else if(*name=="opt_stop")
                        opt_stop=StoppingCondition::from_string(*param);
                    else if(*name=="krylov_stop")
                        krylov_stop=KrylovStop::from_string(*param);
                    else if(*name=="H_type")
                        H_type=Operators::from_string(*param);
                    else if(*name=="Minv_type")
                        Minv_type=Operators::from_string(*param);
                    else if(*name=="dir")
                        dir=LineSearchDirection::from_string(*param);
                    else if(*name=="kind")
                        kind=LineSearchKind::from_string(*param);
                }
            }
            
            // ---------- VERIFYING LABELS ----------- 
        private: 
            // Checks whether we have a valid variable label
            struct is_var : public std::unary_function<std::string, bool> {
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
           
            // Checks whether we have a valid real label
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
                        name == "obj_x" || 
                        name == "obj_xps" ||
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

            // Checks whether we have a valid natural number label
            struct is_nat : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( name == "stored_history" ||
                        name == "history_reset" || 
                        name == "iter" || 
                        name == "iter_max" || 
                        name == "krylov_iter" || 
                        name == "krylov_iter_max" ||
                        name == "krylov_iter_total" || 
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
           
            // Checks whether we have a valid parameter label
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

        protected:
            // Checks that the labels used during serialization are correct
            void checkLabels(
                const Messaging& msg,
                const X_Vectors& xs,
                const Reals& reals,
                const Nats& nats,
                const Params& params
            ) {
                // Create a base message
                const std::string base
                    ="During serialization, found an invalid ";

                // Check the variable names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        xs.first.begin(), xs.first.end(),
                        not1(is_var()));
                    if(name!=xs.first.end()){
                        std::stringstream ss;
                        ss << base << "variable name: " << *name;
                        msg.error(ss.str());
                    }
                }
                
                // Check the real names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        reals.first.begin(), reals.first.end(),
                        not1(is_real()));
                    if(name!=reals.first.end()){
                        std::stringstream ss;
                        ss << base << "real name: " << *name;
                        msg.error(ss.str());
                    }
                }
                    
                // Check the natural names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        nats.first.begin(), nats.first.end(),
                        not1(is_nat()));
                    if(name!=nats.first.end()){
                        std::stringstream ss;
                        ss << base << "natural name: " << *name;
                        msg.error(ss.str());
                    }
                }
                
                // Check the parameter names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        params.first.begin(), params.first.end(),
                        not1(is_param()));
                    if(name!=params.first.end()){
                        std::stringstream ss;
                        ss << base << "parameter name: " << *name;
                        msg.error(ss.str());
                    }
                }
            }

            // Check that the strings used to represent the parameters are
            // correct.
            void checkParams(const Messaging& msg,Params& params) {
                // Create a base message
                const std::string base
                    ="During serialization, found an invalid ";
        
                // Check that the actual parameters are valid
                std::list <std::string>::iterator param=params.second.begin();
                for(std::list <std::string>::iterator name=params.first.begin();
                    name!=params.first.end();
                    name++,param++
                ){
                    // Check the algorithm class
                    if(*name=="algorithm_class"){
                        if(!AlgorithmClass::is_valid()(*param)){
                            std::stringstream ss;
                            ss << base << "algorithm class: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the optimization stopping conditions
                    } else if(*name=="opt_stop"){
                        if(!StoppingCondition::is_valid()(*param)){
                            std::stringstream ss;
                            ss << base << "stopping condition: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the Krylov stopping conditions
                    } else if(*name=="krylov_stop"){
                        if(!KrylovStop::is_valid()(*param)) {
                            std::stringstream ss;
                            ss <<base <<"Krylov stopping condition: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the Hessian type
                    } else if(*name=="H_type"){
                        if(!Operators::is_valid()(*param)) {
                            std::stringstream ss;
                            ss << base<<"Hessian type: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the type of the preconditioner
                    } else if(*name=="Minv_type"){
                        if(!Operators::is_valid()(*param)){
                            std::stringstream ss;
                            ss << base << "preconditioner type: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the line-search direction
                    } else if(*name=="dir"){
                        if(!LineSearchDirection::is_valid()(*param)) {
                            std::stringstream ss;
                            ss << base << "line-search direction: " << *param;
                            msg.error(ss.str());
                        }

                    // Check the kind of line-search
                    } else if(*name=="kind"){
                        if(!LineSearchKind::is_valid()(*param)) {
                            std::stringstream ss;
                            ss << base << "line-search kind: " << *param;
                            msg.error(ss.str());
                        }
                    }
                }
            }

        public: 
            // Release the data into structures controlled by the user 
            void release(
                X_Vectors& xs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::stateToVectors(xs);

                // Copy out all of the scalar information
                Unconstrained <Real,XX>::stateToScalars(reals,nats,params);
            }
            
            // Capture data from structures controlled by the user.  Note,
            // we don't sort the oldY and oldS based on the prefix.  In fact,
            // we completely ignore this information.  Therefore, this routine
            // really depends on oldY and oldS to have their elements inserted
            // into vars in order.  In other words, oldY_1 must come before
            // oldY_2, etc.
            void capture(
                const Messaging& msg,
                X_Vectors& xs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                Unconstrained <Real,XX>::checkLabels(msg,xs,reals,nats,params);

                // Check the strings used to represent parameters
                Unconstrained <Real,XX>::checkParams(msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::vectorsToState(xs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>::scalarsToState(reals,nats,params);

                // Check that we have a valid state 
                Unconstrained <Real,XX>::check(msg);
            }
        };
        
        // State of an equality constrained optimization problem of the form
        // 
        // min_{x \in X} f(x) st g(x) = 0
        //
        // where f : X -> R and g : X -> Y
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY
        > 
        struct EqualityConstrained: public virtual Unconstrained <Real,XX> {
        public:
            // Create some shortcuts for some type names
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;
            
            // The Lagrange multiplier (dual variable) for the equality
            // constraints
            std::list <Y_Vector> y;
            
            // Initialize the state without setting up any variables. 
            EqualityConstrained() {
                // Initialize 
                Unconstrained <Real,XX>::init_params();
                EqualityConstrained <Real,XX,YY>::init_params();
            };
            
            // Initialize the state for equality constrained optimization.
            EqualityConstrained(
                const X_Vector& x,
                const Y_Vector& y
            ) {
                Unconstrained <Real,XX>::init_params();
                Unconstrained <Real,XX>::init_vectors(x);
                EqualityConstrained <Real,XX,YY>::init_params();
                EqualityConstrained <Real,XX,YY>::init_vectors(y);
            }

        protected: 
            // This initializes all the parameters required for equality
            // constrained optimization.  
            void init_params() { }

            // This initializes all the variables required for equality
            // constrained optimization.  
            void init_vectors(const Y_Vector& y_) {
                y.push_back(Y::create()); Y::init(y_,y.back());
                    Y::copy(y_,y.back());
            }

        public:
            // Types for holding information for vectors for restarts
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Y_Vector> > Y_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <unsigned int> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 
        
        protected:
            // Copy out all equality multipliers 
            void stateToVectors(Y_Vectors& ys) {
                ys.first.push_back("y");
                ys.second.splice(ys.second.end(),y);
            }

            // Copy out all the scalar information
            void stateToScalars(
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }
            
            // Copy in all equality multipliers 
            void vectorsToState(Y_Vectors& ys) { 
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
                    if(*name=="y") y.splice(y.end(),ys.second,yy);
                }
            }
            
            // Copy in all the scalar information
            void scalarsToState(
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }

        private:    
            // Checks whether we have a valid equality multiplier label
            struct is_eq : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( name == "y") 
                        return true;
                    else
                        return false;
                }
            };

        protected:
            // Checks that the labels used during serialization are correct
            void checkLabels(
                const Messaging& msg,
                const Y_Vectors& ys,
                const Reals& reals,
                const Nats& nats,
                const Params& params
            ) {
                // Create a base message
                const std::string base
                    ="During serialization, found an invalid ";

                // Check the equality multiplier names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        ys.first.begin(), ys.first.end(),
                        not1(is_eq()));
                    if(name!=ys.first.end()){
                        std::stringstream ss;
                        ss << base << "equality multiplier name: " << *name;
                        msg.error(ss.str());
                    }
                }
            }

            // Check that the strings used to represent the parameters are
            // correct.
            void checkParams(const Messaging& msg,Params& params) { }
            
            // Check that we have a valid set of parameters
            void check(const Messaging& msg){ }

        public: 
            // Release the data into structures controlled by the user 
            void release(
                X_Vectors& xs,
                Y_Vectors& ys,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::stateToVectors(xs);
                EqualityConstrained <Real,XX,YY>::stateToVectors(ys);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>::stateToScalars(reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::stateToScalars(reals,nats,params);
            }

            // Capture data from structures controlled by the user.  
            void capture(
                const Messaging& msg,
                X_Vectors& xs,
                Y_Vectors& ys,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                Unconstrained <Real,XX>::checkLabels(msg,xs,reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::checkLabels(msg,ys,reals,nats,params);

                // Check the strings used to represent parameters
                Unconstrained <Real,XX>::checkParams(msg,params);
                EqualityConstrained <Real,XX,YY>::checkParams(msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::vectorsToState(xs);
                EqualityConstrained <Real,XX,YY>::vectorsToState(ys);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>::scalarsToState(reals,nats,params);
                EqualityConstrained <Real,XX,YY>
                    ::scalarsToState(reals,nats,params);

                // Check that we have a valid state 
                Unconstrained <Real,XX>::check(msg);
                EqualityConstrained <Real,XX,YY>::check(msg);
            }
        };
        
        // State of an inequality constrained optimization problem of the form
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
        public:
            // Create some shortcuts for some type names
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef ZZ <Real> Z;
            typedef typename Z::Vector Z_Vector;
            
            // The Lagrange multiplier (dual variable) for the
            // inequality constraints 
            std::list <Z_Vector> z;
            
            // Initialize the state without setting up any variables. 
            InequalityConstrained() {
                Unconstrained <Real,XX>::init_params();
                InequalityConstrained <Real,XX,ZZ>::init_params();
            };
            
            // Initialize the state for inequality constrained optimization.
            InequalityConstrained(
                const X_Vector& x,
                const Z_Vector& z
            ) {
                Unconstrained <Real,XX>::init_params();
                Unconstrained <Real,XX>::init_vectors(x);
                InequalityConstrained <Real,XX,ZZ>::init_params();
                InequalityConstrained <Real,XX,ZZ>::init_vectors(z);
            }

        protected:
            // This initializes all the parameters required for inequality
            // constrained optimization.  
            void init_params() { }

            // This initializes all the variables required for inequality
            // constrained optimization.  
            void init_vectors(const Z_Vector& z_) {
                z.push_back(Z::create()); Z::init(z_,z.back());
                    Z::copy(z_,z.back());
            }

        public:
            // Types for holding information for vectors for restarts
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Z_Vector> > Z_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <unsigned int> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 

        protected: 
            // Copy out the inequality multipliers 
            void stateToVectors(Z_Vectors& zs) {
                zs.first.push_back("z");
                zs.second.splice(zs.second.end(),z);
            }
            
            // Copy out the scalar information
            void stateToScalars(
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }
            
            // Copy in inequality multipliers 
            void vectorsToState(Z_Vectors& zs) { 
                for(std::list <std::string>::iterator name=zs.first.begin();
                    name!=zs.first.end();
                    name++
                ) {
                    // Since we're using a splice operation, we slowly empty
                    // the multiplier list.  Hence, we always take the first
                    // element.
                    typename std::list <Z_Vector>::iterator zz
                        =zs.second.begin();

                    // Determine which variable we're reading in and then splice
                    // it in the correct location
                    if(*name=="z") z.splice(z.end(),zs.second,zz);
                }
            }
            
            // Copy in the scalar information
            void scalarsToState(
                Reals& reals,
                Nats& nats,
                Params& params
            ) { }


        private:
            // Checks whether we have a valid inequality multiplier label
            struct is_ineq : public std::unary_function<std::string, bool> {
                bool operator () (const std::string& name) const {
                    if( name == "z") 
                        return true;
                    else
                        return false;
                }
            };

        protected:
            // Checks that the labels used during serialization are correct
            void checkLabels(
                const Messaging& msg,
                const Z_Vectors& zs,
                const Reals& reals,
                const Nats& nats,
                const Params& params
            ) {
                // Create a base message
                const std::string base
                    ="During serialization, found an invalid ";

                // Check the inequality multiplier names
                {
                    std::list <std::string>::const_iterator name = find_if(
                        zs.first.begin(), zs.first.end(),
                        not1(is_ineq()));
                    if(name!=zs.first.end()){
                        std::stringstream ss;
                        ss << base << "inequality multiplier name: " << *name;
                        msg.error(ss.str());
                    }
                }
            }

            // Check that the strings used to represent the parameters are
            // correct.
            void checkParams(const Messaging& msg,Params& params) { }
            
            // Check that we have a valid set of parameters.
            void check(const Messaging& msg){ }

        public: 
            // Release the data into structures controlled by the user 
            void release(
                X_Vectors& xs,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,XX>::stateToVectors(xs);
                InequalityConstrained <Real,XX,ZZ>::stateToVectors(zs);
            
                // Copy out all of the scalar information
                Unconstrained <Real,XX>::stateToScalars(reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::stateToScalars(reals,nats,params);
            }
            
            // Capture data from structures controlled by the user.  
            void capture(
                const Messaging& msg,
                X_Vectors& xs,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                Unconstrained <Real,XX>::checkLabels(msg,xs,reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>::checkLabels(
                    msg,zs,reals,nats,params);

                // Check the strings used to represent parameters
                Unconstrained <Real,XX>::checkParams(msg,params);
                InequalityConstrained <Real,XX,ZZ>::checkParams(msg,params);

                // Copy in the variables 
                Unconstrained <Real,XX>::vectorsToState(xs);
                InequalityConstrained <Real,XX,ZZ>::vectorsToState(zs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,XX>::scalarsToState(reals,nats,params);
                InequalityConstrained <Real,XX,ZZ>
                    ::scalarsToState(reals,nats,params);

                // Check that we have a valid state 
                Unconstrained <Real,XX>::check(msg);
                InequalityConstrained <Real,XX,ZZ>::check(msg);
            }
        };

        // State of an equality and inequality constrained optimization
        // problem of the form
        // 
        // min_{x \in X} f(x) st g(x) = 0, h(x) >=_K 0
        //
        // where f : X -> R, g : X -> Y, and h : X -> Z
        template <
            typename Real,
            template <typename> class X,
            template <typename> class Y,
            template <typename> class Z
        > 
        struct Constrained:
            public EqualityConstrained <Real,X,Y>,
            public InequalityConstrained <Real,X,Z>
        {
        public:
            // Create some shortcuts for some type names
            typedef typename X <Real>::Vector X_Vector;
            typedef typename Y <Real>::Vector Y_Vector;
            typedef typename Z <Real>::Vector Z_Vector;
            
            // Initialize the state without setting up any variables. 
            Constrained() {
                Unconstrained <Real,X>::init_params();
                EqualityConstrained <Real,X,Y>::init_params();
                InequalityConstrained <Real,X,Z>::init_params();
            };
            
            // Initialize the state for general constrained optimization.
            Constrained(
                const X_Vector& x,
                const Y_Vector& y,
                const Z_Vector& z
            ) {
                Unconstrained <Real,X>::init_params();
                Unconstrained <Real,X>::init_vectors(x);
                EqualityConstrained <Real,X,Y>::init_params();
                EqualityConstrained <Real,X,Y>::init_vectors(y);
                InequalityConstrained <Real,X,Z>::init_params();
                InequalityConstrained <Real,X,Z>::init_vectors(z);
            }
            
            // Types for holding information for vectors for restarts
            typedef std::pair < std::list <std::string>,
                                std::list <X_Vector> > X_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Y_Vector> > Y_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Z_Vector> > Z_Vectors;
            typedef std::pair < std::list <std::string>,
                                std::list <Real> > Reals;
            typedef std::pair < std::list <std::string>,
                                std::list <unsigned int> > Nats;
            typedef std::pair < std::list <std::string>,
                                std::list <std::string> > Params; 

        public: 
            // Release the data into structures controlled by the user 
            void release(
                X_Vectors& xs,
                Y_Vectors& ys,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {
                // Copy out all of the variable information
                Unconstrained <Real,X>::stateToVectors(xs);
                EqualityConstrained <Real,X,Y>::stateToVectors(ys);
                InequalityConstrained <Real,X,Z>::stateToVectors(zs);
            
                // Copy out all of the scalar information
                Unconstrained <Real,X>::stateToScalars(reals,nats,params);
                EqualityConstrained <Real,X,Y>
                    ::stateToScalars(reals,nats,params);
                InequalityConstrained <Real,X,Z>
                    ::stateToScalars(reals,nats,params);
            }
            
            // Capture data from structures controlled by the user.  
            void capture(
                const Messaging& msg,
                X_Vectors& xs,
                Y_Vectors& ys,
                Z_Vectors& zs,
                Reals& reals,
                Nats& nats,
                Params& params
            ) {

                // Check the labels on the user input
                Unconstrained <Real,X>::checkLabels(msg,xs,reals,nats,params);
                EqualityConstrained <Real,X,Y>::checkLabels(
                    msg,ys,reals,nats,params);
                InequalityConstrained <Real,X,Z>::checkLabels(
                    msg,zs,reals,nats,params);

                // Check the strings used to represent parameters
                Unconstrained <Real,X>::checkParams(msg,params);
                EqualityConstrained <Real,X,Y>::checkParams(msg,params);
                InequalityConstrained <Real,X,Z>::checkParams(msg,params);

                // Copy in the variables 
                Unconstrained <Real,X>::vectorsToState(xs);
                EqualityConstrained <Real,X,Y>::vectorsToState(ys);
                InequalityConstrained <Real,X,Z>::vectorsToState(zs);
                
                // Copy in all of the scalar information
                Unconstrained <Real,X>::scalarsToState(reals,nats,params);
                EqualityConstrained <Real,X,Y>
                    ::scalarsToState(reals,nats,params);
                InequalityConstrained <Real,X,Z>
                    ::scalarsToState(reals,nats,params);

                // Check that we have a valid state 
                Unconstrained <Real,X>::check(msg);
                EqualityConstrained <Real,X,Y>::check(msg);
                InequalityConstrained <Real,X,Z>::check(msg);
            }
        };
    }

    // A function that has free reign to manipulate or analyze the state.
    // This should be used cautiously.
    template <typename State>
    class StateManipulator {
    public:
        // Application
        virtual void operator () (State& state) const {};

        // Allow the derived class to deallocate memory
        virtual ~StateManipulator() {}
    };

    // A simple operator specification, A : X->Y
    template <
        typename Real,
        template <typename> class X,
        template <typename> class Y
    >
    struct Operator {
    private:
        // Create some type shortcuts
        typedef typename X <Real>::Vector X_Vector;
        typedef typename Y <Real>::Vector Y_Vector;
    public:
        // Basic application
        virtual void operator () (const X_Vector& x,Y_Vector &y) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~Operator() {}
    };

    // This contains the different operators required for each class of 
    // optimization methods.
    namespace Operators {

        // Unconstrained optimization 
        template <typename Real,template <typename> class XX>
        struct Unconstrained {
        private:
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector Vector;

        public:

            // The identity operator 
            struct Identity : public Operator <Real,XX,XX> {
                void operator () (const Vector& dx,Vector& result) const{
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
                ScaledIdentity(const State::Unconstrained <Real,XX>& state)
                    : norm_g(state.norm_g), delta_max(state.delta_max) {};

                void operator () (const Vector& dx,Vector& result) const{
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
                const std::list<Vector>& oldY;
                const std::list<Vector>& oldS;

                // Messaging device in case the qusi-Newton information is bad
                const Messaging& msg;
            public:
                BFGS(
                    const Messaging& msg_,
                    const State::Unconstrained <Real,XX>& state
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
                void operator () (const Vector& p, Vector& result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the BFGS Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <Vector> work(oldY.size(),p);
                    for(typename std::list <Vector>::iterator w=work.begin();
                        w!=work.end();
                        w++
                    ) X::init(p,*w);

                    // If we have no vectors in our history, we return the
                    // direction
                    X::copy(p,result);
                    if(oldY.size() == 0) return;

                    // As a safety check, insure that the inner product
                    // between all the (s,y) pairs is positive
                    typename std::list <Vector>::const_iterator y0=oldY.begin();
                    typename std::list <Vector>::const_iterator s0=oldS.begin();
                    while(y0!=oldY.end()){
                        Real inner_y_s=X::innr(*y0++,*s0++);
                        if(inner_y_s<0)
                            msg.error("Detected a (s,y) pair in BFGS that "
                                "possesed a nonpositive inner product");
                    }

                    // Othwerwise, we copy all of the trial step differences
                    // into the work space
                    typename std::list <Vector>::iterator Bisj_iter
                        =work.begin();
                    typename std::list <Vector>::const_iterator sk_iter
                        =oldS.begin();
                    while(Bisj_iter!=work.end())
                        X::copy((*sk_iter++),(*Bisj_iter++));

                    // Keep track of the element Bisi
                    typename std::list <Vector>::const_iterator Bisi_iter
                        =work.end(); Bisi_iter--;

                    // Keep iterating until Bisi equals the first element in the
                    // work list. This means we have computed B1s1, B2s2, ...,
                    // Bksk.
                    Bisj_iter=work.begin();
                    typename std::list <Vector>::const_iterator si_iter
                        =oldS.end(); si_iter--;
                    typename std::list <Vector>::const_iterator yi_iter
                        =oldY.end(); yi_iter--;
                    typename std::list <Vector>::const_iterator sj_iter
                        =oldS.begin();
                    while(1){

                        // Create some reference to our iterators that are
                        // easier to work with
                        const Vector& si=*si_iter;
                        const Vector& yi=*yi_iter;
                        const Vector& Bisi=*Bisi_iter;

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
                            const Vector& sj=*sj_iter;
                            Vector& Bisj=*Bisj_iter;

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
                // Create some type shortcuts
                typedef XX <Real> X;
                typedef typename X::Vector Vector;

                // Stored quasi-Newton information
                const std::list<Vector>& oldY;
                const std::list<Vector>& oldS;

                // Messaging device in case the qusi-Newton information is bad
                const Messaging& msg;
            public:
                SR1(
                    const Messaging& msg_,
                    const State::Unconstrained <Real,XX>& state
                ) : oldY(state.oldY), oldS(state.oldS), msg(msg_) {};
                
                // Operator interface
                void operator () (const Vector& p,Vector& result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the SR1 Hessian approximation, the "
                            "number of stored gradient differences must equal "
                            "the number of stored trial step differences.");

                    // Allocate memory for work
                    std::list <Vector> work(oldY.size(),p);
                    for(typename std::list <Vector>::iterator w=work.begin();
                        w!=work.end();
                        w++
                    ) X::init(p,*w);

                    // If we have no vectors in our history, we return the 
                    // direction
                    X::copy(p,result);
                    if(oldY.size() == 0) return;

                    // Othwerwise, we copy all of the trial step differences 
                    // into the work space
                    typename std::list <Vector>::iterator Bisj_iter
                        =work.begin();
                    typename std::list <Vector>::const_iterator sk_iter
                        =oldS.begin();
                    while(Bisj_iter!=work.end())
                        X::copy((*sk_iter++),(*Bisj_iter++));

                    // Keep track of the element Bisi
                    typename std::list <Vector>::const_iterator Bisi_iter
                        =work.end(); Bisi_iter--;

                    // Keep iterating until Bisi equals the first element in the
                    // work list. This means we have computed B1s1, B2s2, ...,
                    // Bksk.
                    Bisj_iter=work.begin();
                    typename std::list <Vector>::const_iterator si_iter
                        =oldS.end(); si_iter--;
                    typename std::list <Vector>::const_iterator yi_iter
                        =oldY.end(); yi_iter--;
                    typename std::list <Vector>::const_iterator sj_iter
                        =oldS.begin();
                    while(1){

                        // Create some reference to our iterators that are 
                        // easier to work with
                        const Vector& si=*si_iter;
                        const Vector& yi=*yi_iter;
                        const Vector& Bisi=*Bisi_iter;

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
                            const Vector& sj=*sj_iter;
                            Vector& Bisj=*Bisj_iter;

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
                const std::list<Vector>& oldY;
                const std::list<Vector>& oldS;

                // Messaging device in case the qusi-Newton information is bad
                const Messaging& msg;
            public:
                InvBFGS(
                    const Messaging& msg_,
                    const State::Unconstrained <Real,XX>& state
                ) : oldY(state.oldY), oldS(state.oldS), msg(msg_) {};
                
                // Operator interface
                void operator () (const Vector& p,Vector& result) const{

                    // Check that the number of stored gradient and trial step
                    // differences is the same.
                    if(oldY.size() != oldS.size())
                        msg.error("In the inverse BFGS operator, the number "
                            "of stored gradient differences must equal the "
                            "number of stored trial step differences.");
                    
                    // As a safety check, insure that the inner product between
                    // all the (s,y) pairs is positive
                    typename std::list <Vector>::const_iterator y0=oldY.begin();
                    typename std::list <Vector>::const_iterator s0=oldS.begin();
                    while(y0!=oldY.end()){
                        Real inner_y_s=X::innr(*y0++,*s0++);
                        if(inner_y_s<0)
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
                    typename std::list <Vector>::const_iterator y_iter
                        =oldY.begin();
                    typename std::list <Vector>::const_iterator s_iter
                        =oldS.begin();
                    int i=0;
                    while(y_iter != oldY.end()){
                        // Find y_k, s_k, and their inner product
                        const Vector& y_k=*(y_iter++);
                        const Vector& s_k=*(s_iter++);
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
                        const Vector& s_k=*(--s_iter);
                        const Vector& y_k=*(--y_iter);

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
                    const State::Unconstrained <Real,XX>& state
                ) : sr1(msg,state) {};
                void operator () (const Vector& p,Vector& result) const{
                    sr1(p,result);
                }
            };
        };
    }

    // A simple scalar valued function interface, f:X->R
    template <typename Real,template <typename> class XX>
    struct ScalarValuedFunction {
    private:
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector Vector;

        // Hessian approximation
        std::auto_ptr <Operator <Real,XX,XX> > H;

        // This forces derived classes to call the constructor that depends
        // on the state
        ScalarValuedFunction() {}
    public:
         // <- f(x) 
         virtual Real operator () (const Vector& x) const = 0;

         // g = grad f(x) 
         virtual void grad(const Vector& x,Vector& g) const = 0;

         // H_dx = hess f(x) dx 
         virtual void hessvec(const Vector& x,const Vector& dx,Vector& H_dx 
         ) const {
                X::copy(dx,H_dx); 
         }

         // This actually computes the Hessian-vector product.  In essence,
         // we may want to use a Hessian approximation provided by the
         // optimization routines.  The following routine selects whether or
         // not we use the hessvec provided by the user.
         void hv(const Vector& x,const Vector& dx,Vector& H_dx) const {
             if(H.get()!=NULL) 
                (*H)(dx,H_dx);
             else
                hessvec(x,dx,H_dx);
         }

         // The constructor determines whether we really need to build
         // a Hessian-vector product or if we use an internal approximation
         ScalarValuedFunction(
             const Messaging& msg,
             const State::Unconstrained <Real,XX>& state
         ) {
            
            // Determine the Hessian approximation
            switch(state.H_type){
                case Operators::Identity:
                    H.reset(new typename Operators::Unconstrained <Real,XX>
                        ::Identity());
                    break;
                case Operators::ScaledIdentity:
                    H.reset(new typename Operators::Unconstrained <Real,XX>
                        ::ScaledIdentity (state));
                    break;
                case Operators::BFGS:
                    H.reset(new typename Operators::Unconstrained <Real,XX>
                        ::BFGS(msg,state));
                    break;
                case Operators::SR1:
                    H.reset(new typename Operators::Unconstrained <Real,XX>
                        ::SR1(msg,state));
                    break;
                case Operators::External:
                    break;
                default:
                    msg.error("Not a valid Hessian approximation.");
                    break;
            }
        }

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

    // An inequality constraint.  Basically, this is a vector valued function,
    // but it also adds a line-search component.
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class YY 
    >
    struct InequalityConstraint : public VectorValuedFunction <Real,XX,YY> {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector; 
        typedef YY <Real> Y;
        typedef typename Y::Vector Y_Vector; 

        // Line-search, srch <- alpha where max(alpha >=0 : h(x+alpha dx) >=0).
        // In the case where this number is infinite, set alpha=Real(-1.).
        virtual Real srch(const X_Vector& x,const X_Vector& dx) const = 0;

        // We only allow affine functions in inequality constraints.  Hence,
        // the second derivative is zero.
        void pps(
             const X_Vector& x,
             const X_Vector& dx,
             const Y_Vector& dy,
             X_Vector& z
         ) const {
                X::zero(z);
         }
    };

    // All the functions required by an optimization algorithm.  Note, this
    // routine owns the memory for these operations.  
    namespace Functions {
        // Functions pertinent to an unconstrained optimization problem 
        // 
        // min_{x \in X} f(x)
        //
        // where f : X -> R
        template <typename Real,template <typename> class XX> 
        struct Unconstrained {
        public:
            // Objective function
            std::auto_ptr <ScalarValuedFunction <Real,XX> > f;

            // Preconditioner for the Hessian of the objective
            std::auto_ptr <Operator <Real,XX,XX> > Minv;

        protected:
            // Check that all the functions are defined
            void check(const Messaging& msg) {
                // Check that objective function exists 
                if(f.get()==NULL)
                    msg.error("Missing an objective function definition.");
                
                // Check that a preconditioner exists 
                if(Minv.get()==NULL)
                    msg.error("Missing a preconditioner definition.");
            }

            // Create a preconditioner if requested
            void create(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state
            ) {
                // Determine the preconditioner
                switch(state.Minv_type){
                    case Operators::Identity:
                        Minv.reset(new typename
                            Operators::Unconstrained<Real,XX>::Identity());
                        break;
                    case Operators::InvBFGS:
                        Minv.reset(new typename
                            Operators::Unconstrained<Real,XX>
                                ::InvBFGS(msg,state));
                        break;
                    case Operators::InvSR1:
                        Minv.reset(new typename
                            Operators::Unconstrained<Real,XX>
                                ::InvSR1(msg,state));
                        break;
                    case Operators::External:
                        if(Minv.get()==NULL)
                            msg.error("An externally defined preconditioner "
                                "must be provided explicitely.");
                        break;
                    default:
                        msg.error("Not a valid Hessian approximation.");
                        break;
                }
            }

        public:
            // Finalize the construction of the functions by creating
            // any missing pieces and then checking the result.
            void finalize(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state
            ) {

                // Create the missing pieces
                Unconstrained <Real,XX>::create(msg,state);

                // Check the results
                Unconstrained <Real,XX>::check(msg);
            }
        };
        
        // Functions pertinent to an equality constrained optimization problem 
        // 
        // min_{x \in X} f(x) st g(x)=0
        //
        // where f : X -> R, g : X -> Y
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY
        > 
        struct EqualityConstrained: public virtual Unconstrained <Real,XX> {
        public:
            // Equality constraints 
            std::auto_ptr <VectorValuedFunction <Real,XX,YY> > g;
           
        protected:
            // Check that all the functions are defined
            void check(const Messaging& msg) {
                // Check that the constraints exist 
                if(g.get()==NULL)
                    msg.error("Missing the equality constraint definition.");
            }

            // At the moment, we don't create any preconditioners for
            // equality constrained optimization
            void create(
                const Messaging& msg,
                const State::EqualityConstrained <Real,XX,YY>& state
            ) { }

        public:
            // Finalize the construction of the functions by creating
            // any missing pieces and then checking the result.
            void finalize(
                const Messaging& msg,
                const State::EqualityConstrained <Real,XX,YY>& state
            ) {

                // Create the missing pieces
                Unconstrained <Real,XX>::create(msg,state);
                EqualityConstrained <Real,XX,YY>::create(msg,state);

                // Check the results
                Unconstrained <Real,XX>::check(msg);
                EqualityConstrained <Real,XX,YY>::check(msg);
            }
        };
        
        // Functions pertinent to an inequality constrained optimization problem
        // 
        // min_{x \in X} f(x) st h(x) >=_K 0
        //
        // where f : X -> R, h : X -> Z
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class ZZ 
        > 
        struct InequalityConstrained : public virtual Unconstrained <Real,XX> {
        public:
            // Inequality constraints 
            std::auto_ptr <InequalityConstraint <Real,XX,ZZ> > h;

        protected:
            // Check that all the functions are defined
            void check(const Messaging& msg) {
                // Check that the constraints exist 
                if(h.get()==NULL)
                    msg.error("Missing the inequality constraint definition.");
            }

            // At the moment, we don't create any preconditioners for
            // inequality constrained optimization
            void create(
                const Messaging& msg,
                const State::InequalityConstrained <Real,XX,ZZ>& state
            ) { }

        public:
            // Finalize the construction of the functions by creating
            // any missing pieces and then checking the result.
            void finalize(
                const Messaging& msg,
                const State::InequalityConstrained <Real,XX,ZZ>& state
            ) {

                // Create the missing pieces
                Unconstrained <Real,XX>::create(msg,state);
                InequalityConstrained <Real,XX,ZZ>::create(msg,state);

                // Check the results
                Unconstrained <Real,XX>::check(msg);
                InequalityConstrained <Real,XX,ZZ>::check(msg);
            }
        };
        
        // Functions pertinent to a general constrained optimization problem
        // 
        // min_{x \in X} f(x) st g(x)=0, h(x) >=_K 0
        //
        // where f : X -> R, g : X -> Y, h : X -> Z
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY,
            template <typename> class ZZ 
        > 
        struct Constrained:
            public EqualityConstrained <Real,XX,YY>,
            public InequalityConstrained <Real,XX,ZZ>
        {
            // Finalize the construction of the functions by creating
            // any missing pieces and then checking the result.
            void finalize(
                const Messaging& msg,
                const State::Constrained <Real,XX,YY,ZZ>& state
            ) {

                // Create the missing pieces
                Unconstrained <Real,XX>::create(msg,state);
                EqualityConstrained <Real,XX,YY>::create(msg,state);
                InequalityConstrained <Real,XX,ZZ>::create(msg,state);

                // Check the results
                Unconstrained <Real,XX>::check(msg);
                EqualityConstrained <Real,XX,YY>::check(msg);
                InequalityConstrained <Real,XX,ZZ>::check(msg);
            }
        };
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
            typename YY <Real>::Vector& dd
        ){
            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector X_Vector;
            typedef YY <Real> Y;
            typedef typename Y::Vector Y_Vector;

            // Create an element for x+eps dx, x-eps dx, etc. 
            X_Vector x_op_dx; X::init(x,x_op_dx);

            // Create an element for f'(x+eps dx)*dy, etc.
            Y_Vector fps_xopdx_dy; Y::init(dd,fps_xopdx_dy);

            // Zero out the directional derivative
            Y::zero(dd);

            // f'(x+eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            Y::axpy(Real(8.),fps_xopdx_dy,dd);

            // f'(x-eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(-epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            Y::axpy(Real(-8.),fps_xopdx_dy,dd);

            // f'(x+2 eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(Real(2.)*epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            Y::axpy(Real(-1.),fps_xopdx_dy,dd);

            // f'(x-2 eps dx)*dy
            X::copy(x,x_op_dx);
            X::axpy(Real(-2.)*epsilon,dx,x_op_dx);
            f.ps(x_op_dx,dy,fps_xopdx_dy);
            Y::axpy(Real(1.),fps_xopdx_dy,dd);

            // Finish the finite difference calculation 
            Y::scal(Real(1.)/(Real(12.)*epsilon),dd);
        }

        // Performs a finite difference test on the gradient of f where  
        // f : X->R is scalar valued.  In other words, we check grad f using f.
        template <
            typename Real,
            template <typename> class XX
        >
        void gradientCheck(
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
            msg.print("Finite difference test on the gradient.",1);
            for(int i=-2;i<=5;i++){
                Real epsilon=pow(Real(.1),i);
                Real dd=directionalDerivative <> (f,x,dx,epsilon);

                std::stringstream ss;
                if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16)
                    << fabs(dd_grad-dd)/(Real(1e-16)+fabs(dd_grad));
                msg.print(ss.str(),1);
            }
        }
        
        // Performs a finite difference test on the hessian of f where f : X->R
        // is scalar valued.  In other words, we check hess f dx using grad f.
        template <
            typename Real,
            template <typename> class XX
        >
        void hessianCheck(
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
            f.hv(x,dx,hess_f_dx);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the Hessian.",1);
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),hess_f_dx,res);

                // Determine the relative error
                Real rel_err=X::norm(res)/(Real(1e-16)+X::norm(hess_f_dx));

                // Print out the differences
                std::stringstream ss;
                if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str(),1);
            }
        }
        
        // This tests the symmetry of the Hessian.  We accomplish this by
        // comparing <H(x)dx,dxx> to <dx,H(x)dxx>.
        template <
            typename Real,
            template <typename> class XX
        >
        void hessianSymmetryCheck(
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
            f.hv(x,dx,H_x_dx);
            
            // Calculate hess f in the direction dxx.  
            X_Vector H_x_dxx; X::init(x,H_x_dxx);
            f.hv(x,dxx,H_x_dxx);
            
            // Calculate <H(x)dx,dxx>
            Real innr_Hxdx_dxx = X::innr(H_x_dx,dxx);
            
            // Calculate <dx,H(x)dxx>
            Real innr_dx_Hxdxx = X::innr(dx,H_x_dxx);

            // Determine the absolute difference between the two.  This really
            // should be zero.
            Real diff=fabs(innr_Hxdx_dxx-innr_dx_Hxdxx);

            // Send a message with the result
            msg.print("Symmetry test on the Hessian of a scalar valued "
                "function.",1);
            std::stringstream ss;
            ss << "The absolute err. between <H(x)dx,dxx> and <dx,H(x)dxx>: "
                << std::scientific << std::setprecision(16) << diff;
            msg.print(ss.str(),1);
        }

        // Performs a finite difference test on the derivative of a
        // vector-valued function f.  Specifically, we check f'(x)dx using f.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        void derivativeCheck(
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
                "vector-valued function.",1);
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,epsilon,res);

                // Determine the residual.  Store in res.
                Y::axpy(Real(-1.),fp_x_dx,res);

                // Determine the relative error
                Real rel_err=Y::norm(res)/(Real(1e-16)+Y::norm(fp_x_dx));

                // Print out the differences
                std::stringstream ss;
                if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str(),1);
            }
        }

        // Performs an adjoint check on the first-order derivative of a vector
        // valued function.  In other words, we check that
        // <f'(x)dx,dy> = <dx,f'(x)*dy>
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        void derivativeAdjointCheck(
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
            msg.print("Adjoint test on the first derivative of a vector valued "
                "function.",1);
            std::stringstream ss;
            ss << "The absolute err. between <f'(x)dx,dy> and <dx,f'(x)*dy>: "
                << std::scientific << std::setprecision(16) << diff;
            msg.print(ss.str(),1);
        }

        // Performs a finite difference test on the second-derivative-adjoint 
        // of a vector-valued function f.  Specifically, we check
        // (f''(x)dx)*dy using f'(x)*dy.
        template <
            typename Real,
            template <typename> class XX,
            template <typename> class YY 
        >
        void secondDerivativeCheck(
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
            X_Vector fpps_x_dx_dy; X::init(dy,fpps_x_dx_dy);
            f.pps(x,dx,dy,fpps_x_dx_dy);

            // Compute an ensemble of finite difference tests in a linear manner
            msg.print("Finite difference test on the 2nd-derivative adj. "
                "of a vector-valued function.",1);
            for(int i=-2;i<=5;i++){

                // Calculate the directional derivative
                Real epsilon=pow(Real(.1),i);
                directionalDerivative <> (f,x,dx,dy,epsilon,res);

                // Determine the residual.  Store in res.
                X::axpy(Real(-1.),fpps_x_dx_dy,res);

                // Determine the relative error
                Real rel_err=X::norm(res)/(Real(1e-16)+X::norm(fpps_x_dx_dy));

                // Print out the differences
                std::stringstream ss;
                if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
                else ss << "The relative difference (1e-" << i << "): ";
                ss << std::scientific << std::setprecision(16) << rel_err; 
                msg.print(ss.str(),1);
            }
        }
    }

    // This contains the different algorithms used for optimization 
    namespace Algorithms {

        // Unconstrained optimization 
        template <
            typename Real,
            template <typename> class XX
        >
        struct Unconstrained {

            // Create some type shortcuts
            typedef XX <Real> X;
            typedef typename X::Vector Vector;

            // Checks a set of stopping conditions
            static StoppingCondition::t checkStop(
                const State::Unconstrained <Real,XX>& state
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

                // Check whether the norm is small relative to some typical
                // gradient
                if(norm_g < eps_g*norm_gtyp)
                    return StoppingCondition::RelativeGradientSmall;

                // Check whether the change in the step length has become too
                // small relative to some typical step
                if(norm_s < eps_s*norm_styp)
                    return StoppingCondition::RelativeStepSmall;

                // Check if we've exceeded the number of iterations
                if(iter>=iter_max)
                    return StoppingCondition::MaxItersExceeded;

                // Otherwise, return that we're not converged 
                return StoppingCondition::NotConverged;
            }
        
            // Computes the truncated-CG (Steihaug-Toint) trial step for
            // trust-region and line-search algorithms
            static void truncatedCG(
                const Messaging& msg,
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){

                // Create shortcuts to some elements in the state
                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const Vector& x=*(state.x.begin());
                const Vector& g=*(state.g.begin());
                const Real& delta=state.delta;
                const Real& eps_cg=state.eps_krylov;
                const unsigned int& iter_max=state.krylov_iter_max;
                Vector& s_k=*(state.s.begin());
                unsigned int& iter=state.krylov_iter;
                unsigned int& iter_total=state.krylov_iter_total;
                KrylovStop::t& krylov_stop=state.krylov_stop;
                Real& rel_err=state.krylov_rel_err;

                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);
                const Operator <Real,XX,XX>& Minv=*(fns.Minv);

                // Allocate memory for temporaries that we need
                Vector g_k; X::init(x,g_k);
                Vector v_k; X::init(x,v_k);
                Vector p_k; X::init(x,p_k);
                Vector H_pk; X::init(x,H_pk);

                // Allocate memory for a few constants that we need to track 
                Real kappa;
                Real sigma;
                Real alpha(0.);
                Real beta;
                Real norm_sk_M2,norm_skp1_M2(0.),norm_pk_M2,norm_g;
                Real inner_sk_M_pk,inner_gk_vk,inner_gkp1_vkp1;

                // Initialize our variables
                X::scal(Real(0.),s_k);       // s_0=0
                X::copy(g,g_k);              // g_0=g
                Minv(g_k,v_k);               // v_0=inv(M)*g_0
                X::copy(v_k,p_k);            // p_0=-v_0
                X::scal(Real(-1.),p_k);
                norm_sk_M2=Real(0.);         // || s_0 ||_M^2 = 0
                norm_pk_M2=X::innr(g_k,v_k); // || p_0 ||_M^2 = <g_0,v_0>
                inner_sk_M_pk=Real(0.);      // <s_0,M p_0>=0
                inner_gk_vk=norm_pk_M2;      // <g_0,v_0> = || p_0 ||_M^2
                norm_g=X::innr(g,g);         // || g ||

                // Run truncated CG until we hit our max iteration or we
                // converge
                iter_total++;
                for(iter=1;iter<=iter_max;iter++,iter_total++){
                    // H_pk=H p_k
                    f.hv(x,p_k,H_pk);

                    // Compute the curvature for this direction.
                    // kappa=<p_k,H p_k>
                    kappa=X::innr(p_k,H_pk);

                    // If we have negative curvature, don't bother with the
                    // next two steps since we're going to exit and we won't
                    // need them.  
                    if(kappa > 0){
                        // Determine a trial point
                        alpha = X::innr(g_k,v_k)/kappa;

                        // || s_k+alpha_k p_k ||
                        norm_skp1_M2=norm_sk_M2+Real(2.)*alpha*inner_sk_M_pk
                            +alpha*alpha*norm_pk_M2;
                    }

                    // If we have negative curvature or our trial point is
                    // outside the trust region radius, terminate truncated-CG
                    // and find our final step.  We ignore the trust-region
                    // radius in the case of a line-search algorithm. We have
                    // the kappa!=kappa check in order to trap NaNs.
                    if( kappa <= 0 ||
                        (algorithm_class==AlgorithmClass::TrustRegion
                            && norm_skp1_M2 >= delta*delta) ||
                        kappa!=kappa
                    ){
                        switch(algorithm_class){
                        // In the case of a trust-region algorithm, take the
                        // last trial step and extend it all of the way to the
                        // trust-region radius.
                        case AlgorithmClass::TrustRegion:
                            // sigma=positive root of
                            // || s_k + sigma p_k ||_M=delta
                            sigma=(-inner_sk_M_pk
                                + sqrt(inner_sk_M_pk*inner_sk_M_pk
                                + norm_pk_M2*(delta*delta-norm_sk_M2)))
                                / norm_pk_M2;

                            // s_kp1=s_k+sigma p_k
                            X::axpy(sigma,p_k,s_k);

                            // Update the residual error for out output,
                            // g_k=g_k+sigma Hp_k
                            X::axpy(sigma,H_pk,g_k);
                            break;
                            
                        // In the case of a line-search algorithm, do nothing
                        // unless we're on the first iteration.  In this case,
                        // take the steepest descent direction. 
                        case AlgorithmClass::LineSearch:
                            if(iter==1){
                                X::copy(g,s_k);
                                X::scal(Real(-1.),s_k);
                            }
                            break;
                        }

                        // Return a message as to why we exited
                        if(kappa<=0 || kappa!=kappa)
                            krylov_stop = KrylovStop::NegativeCurvature;
                        else
                            krylov_stop = KrylovStop::TrustRegionViolated;

                        // Exit the loop
                        break;
                    }

                    // Take a step in the computed direction. s_k=s_k+alpha p_k
                    X::axpy(alpha,p_k,s_k);

                    // Update the norm of sk
                    norm_sk_M2=norm_skp1_M2;
                    
                    // g_k=g_k+alpha H p_k
                    X::axpy(alpha,H_pk,g_k);

                    // Test whether we've converged CG
                    rel_err=sqrt(X::innr(g_k,g_k)) / (Real(1e-16)+norm_g);
                    if(rel_err <= eps_cg){
                        krylov_stop = KrylovStop::RelativeErrorSmall;
                        break;
                    }

                    // v_k = Minv g_k
                    Minv(g_k,v_k);

                    // Compute the new <g_kp1,v_kp1>
                    inner_gkp1_vkp1=X::innr(g_k,v_k);

                    // beta = <g_kp1,v_kp1> / <g_k,v_k>
                    beta= inner_gkp1_vkp1 / inner_gk_vk;

                    // Store the new inner product between g_k and p_k
                    inner_gk_vk=inner_gkp1_vkp1;
                    
                    // Find the new search direction.  p_k=-v_k + beta p_k
                    X::scal(beta,p_k);
                    X::axpy(Real(-1.),v_k,p_k);

                    // Update the inner product between s_k and M p_k
                    inner_sk_M_pk=beta*(inner_sk_M_pk+alpha*norm_pk_M2);

                    // Update the norm of p_k
                    norm_pk_M2=inner_gk_vk+beta*beta*norm_pk_M2; 

                    // Print out diagnostics
                    printKrylov(msg,state);
                }

                // Check if we've exceeded the maximum iteration
                if(iter>iter_max){
                    krylov_stop=KrylovStop::MaxItersExceeded;
                    iter--; iter_total--;
                }
               
                // Grab the relative error in the CG solution.
                // If we're running a line-search algorithm and we exit on the
                // first iteration, then this doesn't make sense since we
                // took the steepest descent direction without the
                // preconditioner.  Hence, we return a NaN.
                if(algorithm_class==AlgorithmClass::LineSearch && iter==1)
                    rel_err= Real(std::numeric_limits<double>::quiet_NaN()); 

                // Otherwise, we calculate the current relative error.
                else
                    rel_err=sqrt(X::innr(g_k,g_k)) / (Real(1e-16)+norm_g);
                    
                // Print out diagnostics
                if(iter!=iter_max) printKrylov(msg,state);
            }

            // Checks whether we accept or reject a step
            static bool checkStep(
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create shortcuts to some elements in the state
                const Vector& s=*(state.s.begin());
                const Vector& x=*(state.x.begin());
                const Vector& g=*(state.g.begin());
                const Real& eta1=state.eta1;
                const Real& eta2=state.eta2;
                const Real& delta_max=state.delta_max;
                const Real& obj_x=state.obj_x;
                const Real& norm_s=state.norm_s;
                Real& delta=state.delta;
                Real& rho=state.rho;
                Real& obj_xps=state.obj_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);

                // Allocate memory for temporaries that we need
                Vector xps; X::init(x,xps);
                Vector Hx_s; X::init(x,Hx_s);

                // Determine x+s 
                X::copy(s,xps);
                X::axpy(Real(1.),x,xps);

                // Determine the objective function evaluated at u+s
                obj_xps=f(xps);
                
                // Determine H(x)s
                f.hv(x,s,Hx_s);

                // Determine alpha+<g,s>+.5*<H(u)s,s>
                Real model_s=obj_x+X::innr(g,s)+Real(.5)*X::innr(Hx_s,s);

                // Add a safety check in case we don't actually minimize the TR
                // subproblem correctly. This could happen for a variety of
                // reasons.  Most notably, if we do not correctly calculate the
                // Hessian approximation, we could have a nonsymmetric 
                // approximation.  In that case, truncated-CG will exit, but 
                // has an undefined result.  In the case that the actual 
                // reduction also increases, rho could have an extraneous 
                // positive value.  Hence, we require an extra check.
                if(model_s > obj_x){
                    delta = norm_s/Real(2.);
                    rho = Real(std::numeric_limits<double>::quiet_NaN()); 
                    return false;
                }

                // Determine the ratio of reductions
                rho = (obj_x - obj_xps) / (obj_x - model_s);

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
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create some shortcuts
                unsigned int& rejected_trustregion=state.rejected_trustregion;
                Vector& s=*(state.s.begin());
                Real& norm_s=state.norm_s;
                std::list <Vector>& oldY=state.oldY; 
                std::list <Vector>& oldS=state.oldS; 
                unsigned int& history_reset=state.history_reset;

                // Continue to look for a step until one comes back as valid
                for(rejected_trustregion=0;
                    true; 
                    rejected_trustregion++
                ) {
                    // If the number of rejected steps is above the
                    // history_reset threshold, destroy the quasi-Newton
                    // information
                    if(rejected_trustregion > history_reset){
                        oldY.clear();
                        oldS.clear();
                    }

                    // Print out diagnostic information if we reject a step
                    if(rejected_trustregion>0) printState(msg,state,true);

                    // Use truncated-CG to find a new trial step
                    truncatedCG(msg,fns,state);

                    // Save the length of the trial step
                    norm_s=sqrt(X::innr(s,s));

                    // Check whether the step is good
                    if(checkStep(fns,state)) break;
                } 
            }
        
            // Steepest descent search direction
            static void SteepestDescent(
                State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                Vector& s=*(state.s.begin());

                // We take the steepest descent direction
                X::copy(g,s);
                X::scal(Real(-1.),s);
            }
    
            // Nonlinear Conjugate Gradient
            static void NonlinearCG(
                const NonlinearCGDirections::t dir,
                State::Unconstrained <Real,XX>& state
            ) {
            
                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                const Vector& s_old=*(state.s_old.begin());
                const int& iter=state.iter;
                Vector& s=*(state.s.begin());

                // If we're on the first iterations, we take the steepest
                // descent direction
                if(iter==1) SteepestDescent(state);

                // On subsequent iterations, we take the specified direction
                else {
                    // Find the momentum parameter
                    double beta;
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
                State::Unconstrained <Real,XX>& state
            ) {

                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                const Vector& g_old=*(state.g_old.begin());

                // Return the momentum parameter
                return X::innr(g,g)/X::innr(g_old,g_old);
            }
        
            // Polak-Ribiere CG search direction
            static Real PolakRibiere(
                State::Unconstrained <Real,XX>& state
            ) {

                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                const Vector& g_old=*(state.g_old.begin());
                    
                // Return the momentum parameter
                return (X::innr(g,g)-X::innr(g,g_old))
                    /X::innr(g_old,g_old);
            }
            
            // Hestenes-Stiefel search direction
            static Real HestenesStiefel(
                State::Unconstrained <Real,XX>& state
            ) {

                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                const Vector& g_old=*(state.g_old.begin());
                const Vector& s_old=*(state.s_old.begin());
                    
                // Return the momentum parameter
                return (X::innr(g,g)-X::innr(g,g_old))
                    /(X::innr(g,s_old)-X::innr(g_old,s_old));
            }

            // BFGS search direction
            static void BFGS(
                const Messaging& msg,
                State::Unconstrained <Real,XX>& state
            ) {
                
                // Create some shortcuts 
                const Vector& g=*(state.g.begin());
                Vector& s=*(state.s.begin());

                // Create the inverse BFGS operator
                typename Operators::Unconstrained <Real,XX>::InvBFGS
                    Hinv(msg,state); 

                // Apply the inverse BFGS operator to the gradient
                Hinv(g,s);

                // Negate the result
                X::scal(Real(-1.),s);
            }

            // Compute a Golden-Section search between eps and 2*alpha where
            // alpha is the last line search parameter.
            static void goldenSection(
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts
                const Vector& x=*(state.x.begin());
                const unsigned int& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& obj_xps=state.obj_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);

                // Create one work element that holds x+mu s or x+lambda s
                Vector x_p_s; X::init(x,x_p_s);

                // Find 1 over the golden ratio
                Real beta=Real(2./(1.+sqrt(5.)));

                // Find a bracket for the linesearch such that a < b
                Real a=Real(1e-16);
                Real b=Real(2.)*alpha;

                // Find two new points between a and b, mu and lambda,
                // such that lambda < mu
                double lambda=a+(1.-beta)*(b-a);
                double mu=a+beta*(b-a);

                // Find the objective value at mu and labmda 

                // mu 
                X::copy(x,x_p_s);
                X::axpy(mu,s,x_p_s);
                Real obj_mu=f(x_p_s);

                // lambda
                X::copy(x,x_p_s);
                X::axpy(lambda,s,x_p_s);
                Real obj_lambda=f(x_p_s);

                // Search for a fixed number of iterations 
                for(iter=0;iter<iter_max;iter++){

                    // If the objective is greater on the left, bracket on the
                    // right.  Alternatively, it's possible that we're going to
                    // generate a NaN on the right.  This means that obj_mu=NaN.
                    // In this case we want to bracket on the left.  Since
                    // obj_lambda > obj_mu will return false when obj_mu is a
                    // NaN, we should be safe.
                    if(obj_lambda > obj_mu){
                        a=lambda;
                        lambda=mu;
                        obj_lambda=obj_mu;
                        mu=a+beta*(b-a);

                        X::copy(x,x_p_s);
                        X::axpy(mu,s,x_p_s);
                        obj_mu=f(x_p_s);

                    // Otherwise, the objective is greater on the right, so
                    // bracket on the left
                    } else {
                        b=mu;
                        mu=lambda;
                        obj_mu=obj_lambda;
                        lambda=a+(1-beta)*(b-a);
                
                        X::copy(x,x_p_s);
                        X::axpy(lambda,s,x_p_s);
                        obj_lambda=f(x_p_s);
                    }
                }

                // Keep track of the total number of iterations
                iter_total += iter;

                // Once we're finished narrowing in on a solution, take our best
                // guess for the line search parameter
                alpha=obj_lambda < obj_mu ? lambda : mu;

                // Save the objective value at this step
                obj_xps=obj_lambda < obj_mu ? obj_lambda : obj_mu;
            }

            // Find the line search parameter based on the 2-point approximation
            // from Barzilai and Borwein
            static void twoPoint(
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts
                const Vector& x=*(state.x.begin());
                const Vector& g=*(state.g.begin());
                const Vector& x_old=*(state.x_old.begin());
                const Vector& g_old=*(state.g_old.begin());
                const LineSearchKind::t& kind=state.kind;
                Real& alpha=state.alpha;
                Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& obj_xps=state.obj_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);

                // Create elements for delta_x and delta_g as well as one work
                // element for storing x+alpha s
                Vector delta_x; X::init(x,delta_x);
                Vector delta_g; X::init(x,delta_g);
                Vector x_p_s; X::init(x,x_p_s);

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

                // Save the objective value at this step
                X::copy(x,x_p_s);
                X::axpy(alpha,s,x_p_s);
                obj_xps=f(x_p_s);

                // Since we do one function evaluation, increase the linesearch
                // iteration by one
                iter=1; iter_total++;
            }
            
            // Compute a backtracking line-search. 
            static void backTracking(
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts
                const Vector& x=*(state.x.begin());
                const unsigned int& iter_max=state.linesearch_iter_max;
                Real& alpha=state.alpha;
                Vector& s=*(state.s.begin());
                unsigned int& iter_total=state.linesearch_iter_total;
                unsigned int& iter=state.linesearch_iter;
                Real& obj_xps=state.obj_xps;
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);

                // Create one work element for holding x+alpha s
                Vector x_p_s; X::init(x,x_p_s);

                // Store the best objective value and alpha that we used to
                // find it.
                // Our initial guess will be at alpha*2.
                Real alpha_best=Real(2.)*alpha;
                X::copy(x,x_p_s);
                X::axpy(alpha_best,s,x_p_s);
                Real obj_best=f(x_p_s);

                // Evaluate the objective iter_max times at a distance of
                // 2*alpha, alpha, alpha/2, ....  Then, pick the best one.
                // Note, we start iter at 1 since we've already done one
                // iteration above.
                Real alpha0=alpha;
                for(iter=1;iter<iter_max;iter++){
                    // Evaluate f(x+alpha*s)
                    X::copy(x,x_p_s);
                    X::axpy(alpha0,s,x_p_s);
                    Real obj=f(x_p_s);

                    // If this is better than our best guess so far, save it
                    if(obj<obj_best){
                        obj_best=obj;
                        alpha_best=alpha0;
                    }

                    // Reduce the size of alpha
                    alpha0 /= Real(2.);
                }

                // Save the best objective and alpha found
                alpha=alpha_best;
                obj_xps=obj_best;

                // Indicate how many iterations we used to find this value
                iter_total+=iter;
            }

            // Finds a trial step using a line-search for globalization
            static void getStepLS(
                const Messaging& msg,
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create some shortcuts
                const LineSearchDirection::t& dir=state.dir;
                const LineSearchKind::t& kind=state.kind;
                const int& iter=state.iter;
                const int& linesearch_iter_max=state.linesearch_iter_max;
                const Real& obj_x=state.obj_x;
                Real& obj_xps=state.obj_xps;
                Vector& s=*(state.s.begin());
                Real& norm_s=state.norm_s;
                Real& alpha=state.alpha;

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
                case LineSearchDirection::NewtonCG:
                    truncatedCG(msg,fns,state);
                    break;
                }

                // Do a line-search in the specified direction
                switch(kind){
                case LineSearchKind::GoldenSection:
                    // Continue doing a line-search until we get a reduction
                    // in the objective value.
                    do {
                        // Conduct the golden section search
                        goldenSection(fns,state);

                        // If we have no reduction in the objective, print
                        // some diagnostic information.
                        if(obj_xps > obj_x) {
                            norm_s=alpha*X::norm(s);
                            printState(msg,state,true);

                            // We reduce alpha by a factor of four when we
                            // reject the step since the line-search always
                            // looks twice alpha out in the direction of the 
                            // search direction.  By reducing alpha by a factor
                            // of four we insure that the next line-search
                            // examines a unique set of points.
                            alpha /= Real(4.);
                        }

                    // If we don't decrease the objective, try again 
                    } while(obj_x < obj_xps || obj_xps!=obj_xps);
                    break;
                case LineSearchKind::BackTracking:
                    // Continue doing a line-search until we get a reduction
                    // in the objective value.
                    do {
                        // Conduct the golden section search
                        goldenSection(fns,state);

                        // If we have no reduction in the objective, print
                        // some diagnostic information.
                        if(obj_xps > obj_x) {
                            norm_s=alpha*X::norm(s);
                            printState(msg,state,true);

                            // We set alpha to be four times less than the
                            // minimimum alpha we searched before.  We do this
                            // since the line-search always looks twice alpha
                            // out in the beginning of the search.
                            alpha = alpha/pow(Real(2.),linesearch_iter_max+1);
                        }

                    // If we don't decrease the objective, try again 
                    } while(obj_x < obj_xps || obj_xps!=obj_xps);
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
                const Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;

                // Choose whether we use a line-search or trust-region method
                switch(algorithm_class){
                case AlgorithmClass::TrustRegion:
                    getStepTR(msg,fns,state);
                    break;
                case AlgorithmClass::LineSearch:
                    getStepLS(msg,fns,state);
                    break;
                }
            }

            // Updates the quasi-Newton information
            static void updateQuasi(
                State::Unconstrained <Real,XX>& state
            ){
                // Exit immediately if we're not using a quasi-Newton method
                if(state.stored_history==0) return;

                // Create some shortcuts
                const Vector& x=*(state.x.begin());
                const Vector& g=*(state.g.begin());
                const Vector& x_old=*(state.x_old.begin());
                const Vector& g_old=*(state.g_old.begin());
                const Operators::t& Minv_type=state.Minv_type;
                const Operators::t& H_type=state.H_type;
                const LineSearchDirection::t& dir=state.dir;
                std::list <Vector>& oldY=state.oldY;
                std::list <Vector>& oldS=state.oldS;
               
                // Allocate some temp storage for y and s
                Vector s; X::init(x,s);
                Vector y; X::init(x,y);

                // Find s = u-u_old
                X::copy(x,s);
                X::axpy(Real(-1.),x_old,s);

                // Find y = g - g_old
                X::copy(g,y);
                X::axpy(Real(-1.),g_old,y);

                // If we're using BFGS, check that <y,s> > 0
                if((Minv_type==Operators::InvBFGS || H_type==Operators::BFGS
                    || dir==LineSearchDirection::BFGS)
                    && X::innr(y,s) < 0)
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
            static void getMin(
                const Messaging& msg,
                const StateManipulator<State::Unconstrained <Real,XX> >& smanip,
                Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create some shortcuts
                Vector& x=*(state.x.begin());
                Vector& g=*(state.g.begin());
                Vector& s=*(state.s.begin());
                Vector& x_old=*(state.x_old.begin());
                Vector& g_old=*(state.g_old.begin());
                Vector& s_old=*(state.s_old.begin());
                Real& obj_x=state.obj_x;
                Real& obj_xps=state.obj_xps;
                Real& norm_s=state.norm_s;
                Real& norm_g=state.norm_g;
                Real& norm_gtyp=state.norm_gtyp;
                Real& norm_styp=state.norm_styp;
                unsigned int& iter=state.iter;
                StoppingCondition::t& opt_stop=state.opt_stop;
                AlgorithmClass::t& algorithm_class=state.algorithm_class;
                LineSearchDirection::t& dir=state.dir;

                // Finalize the functions that define the problem 
                fns.finalize(msg,state);
                
                // Create shortcuts to the functions that we need
                const ScalarValuedFunction <Real,XX>& f=*(fns.f);

                // Evaluate the objective function and gradient if we've not
                // done so already
                if(obj_x != obj_x) {
                    obj_x=f(x);
                    f.grad(x,g);
                    norm_g=sqrt(X::innr(g,g));
                    norm_gtyp=norm_g;
                }

                // Prints out the header for the diagonstic information
                printStateHeader(msg,state);
                if(algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG)
                    printKrylovHeader(msg,state);

                // Primary optimization loop
                do{
                    // Print some diagnostic information
                    printState(msg,state);

                    // Get a new optimization iterate.  
                    getStep(msg,fns,state);

                    // If we've not calculated it already, save the size of
                    // the step
                    if(norm_styp!=norm_styp) norm_styp=norm_s;

                    // Save the old variable, gradient, and trial step.  This
                    // is useful for both CG and quasi-Newton methods.
                    X::copy(x,x_old);
                    X::copy(g,g_old);
                    X::copy(s,s_old);

                    // Move to the new iterate
                    X::axpy(Real(1.),s,x);

                    // Manipulate the state if required
                    smanip(state);

                    // Find the new objective function and gradient
                    obj_x=obj_xps;
                    f.grad(x,g);
                    norm_g=sqrt(X::innr(g,g));

                    // Update the quasi-Newton information
                    updateQuasi(state);

                    // Increase the iteration
                    iter++;

                    // Check the stopping condition
                    opt_stop=checkStop(state);
                } while(opt_stop==StoppingCondition::NotConverged);
                        
                // Print a final diagnostic 
                printState(msg,state);
            }
            
            // Solves an optimization problem where the user doesn't know about
            // the state manipulator
            static void getMin(
                const Messaging& msg,
                Functions::Unconstrained <Real,XX>& fns,
                State::Unconstrained <Real,XX>& state
            ){
                // Create an empty state manipulator
                StateManipulator<State::Unconstrained <Real,XX> > smanip;

                // Minimize the problem
                getMin(msg,smanip,fns,state);
            }
        
            // Prints out useful information regarding the current optimization
            // state
            static void printState(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state,
                const bool noiter=false
            ) {
                // Create some shortcuts
                const int& iter=state.iter;
                const Real& obj_x=state.obj_x;
                const Real& norm_g=state.norm_g;
                const Real& norm_s=state.norm_s;
                const Real& krylov_rel_err=state.krylov_rel_err;
                const int& krylov_iter=state.krylov_iter;
                const KrylovStop::t& krylov_stop=state.krylov_stop; 
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;
                const int& linesearch_iter=state.linesearch_iter;

                // Basic information
                std::stringstream ss;
                if(!noiter) ss << std::setw(4) << iter << ' ';
                else ss << std::setw(4) << '*' << ' ';
                ss << std::scientific << std::setprecision(3)
                    << std::setw(11) << obj_x << ' '
                    << std::setw(11) << norm_g << ' ';
                if(iter==0) ss << "            ";
                else ss << std::setw(11) << norm_s << ' ';

                // Information for the Krylov method
                if(algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                    ss << std::setw(11) << krylov_rel_err << ' ' 
                        << std::setw(6) << krylov_iter << ' ';

                    switch(krylov_stop){
                    case KrylovStop::NegativeCurvature:
                        ss << std::setw(10) << "Neg Curv" << ' '; 
                        break;
                    case KrylovStop::RelativeErrorSmall:
                        ss << std::setw(10) << "Rel Err " << ' '; 
                        break;
                    case KrylovStop::MaxItersExceeded:
                        ss << std::setw(10) << "Max Iter" << ' '; 
                        break;
                    case KrylovStop::TrustRegionViolated:
                        ss << std::setw(10) << "Trst Reg" << ' '; 
                        break;
                    }
                }

                // Information for the line search
                if(algorithm_class==AlgorithmClass::LineSearch){
                    ss << std::setw(6) << linesearch_iter << ' ';
                }

                // Send the information to the screen
                msg.print(ss.str(),1);
            }
            
            // Prints out useful information regarding the Krylov method 
            static void printKrylov(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts
                const Real& krylov_rel_err=state.krylov_rel_err;
                const int& krylov_iter=state.krylov_iter;
                const int& krylov_iter_total=state.krylov_iter_total;

                // Basic information
                std::stringstream ss;
                ss << "  " << std::setw(4) << krylov_iter << ' '
                    << std::setw(6) << krylov_iter_total << ' '
                    << std::scientific << std::setprecision(3)
                        << std::setw(11) << krylov_rel_err; 

                // Send the information to the screen
                msg.print(ss.str(),2);
            }

            // Prints out a description for the state
            static void printStateHeader(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state
            ) {
                // Create some shortcuts
                const AlgorithmClass::t& algorithm_class=state.algorithm_class;
                const LineSearchDirection::t& dir=state.dir;

                // Basic Header
                std::stringstream ss;
                ss << "Iter" << ' '
                        << std::setw(11) << "Obj Value" << ' '
                        << std::setw(11) << "  norm(g)" << ' '
                        << std::setw(11) << "  norm(s)" << ' ';

                // In case we're using a Krylov method
                if(algorithm_class==AlgorithmClass::TrustRegion
                    || dir==LineSearchDirection::NewtonCG
                ){
                        ss << std::setw(11) << "Kry Error" << ' '
                        << std::setw(6) << "KryIt" << ' ' 
                        << std::setw(10) << "Kry Why " << ' ';
                }

                // Information for the line search
                if(algorithm_class==AlgorithmClass::LineSearch){
                    ss << std::setw(6) << "LS It" << ' ';
                }

                // Send the information to the screen
                msg.print(ss.str(),1);
            }
            
            // Prints out a description for the Krylov method 
            static void printKrylovHeader(
                const Messaging& msg,
                const State::Unconstrained <Real,XX>& state
            ) {
                // Basic Header
                std::stringstream ss;
                ss << "  Iter" << ' '
                        << std::setw(6) << "Total" << ' '
                        << std::setw(11) << "Rel Err" << ' ';

                // Send the information to the screen
                msg.print(ss.str(),2);
            }
        
        };
    }
}
#endif
