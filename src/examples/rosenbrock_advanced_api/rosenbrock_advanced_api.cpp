// In this example, we duplicate the Rosenbrock example while demonstrating
// some of the more advanced API features such as custom vector spaces,
// messaging objects, and restarts.

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include "optizelle/optizelle.h"
#include "optizelle/json.h"

// Grab Optizelle's Natural type
using Optizelle::Natural;

//---VectorSpace0---
// Defines the vector space used for optimization.
template <typename Real>
struct MyVS { 
    typedef std::vector <Real> Vector;
    
    // Memory allocation and size setting
    static Vector init(Vector const & x) {
        return std::move(Vector(x.size()));
    }

    // y <- x (Shallow.  No memory allocation.)
    static void copy(Vector const & x, Vector & y) {
        for(Natural i=0;i<x.size();i++){
            y[i]=x[i];
        }
    }

    // x <- alpha * x
    static void scal(const Real& alpha, Vector & x) {
        for(Natural i=0;i<x.size();i++){
            x[i]=alpha*x[i];
        }
    }

    // x <- 0 
    static void zero(Vector & x) {
        for(Natural i=0;i<x.size();i++){
            x[i]=0.;
        }
    }

    // y <- alpha * x + y
    static void axpy(const Real& alpha, Vector const & x, Vector & y) {
        for(Natural i=0;i<x.size();i++){
            y[i]=alpha*x[i]+y[i];
        }
    }

    // innr <- <x,y>
    static Real innr(Vector const & x,Vector const & y) {
        Real z=0;
        for(Natural i=0;i<x.size();i++)
            z+=x[i]*y[i];
        return z;
    }

    // x <- random
    static void rand(Vector & x){
        std::mt19937 gen(1);
        std::uniform_real_distribution<Real> dis(Real(0.),Real(1.));
        for(Natural i=0;i<x.size();i++) 
            x[i]=Real(dis(gen));
    }
    // Jordan product, z <- x o y.
    static void prod(Vector const & x, Vector const & y, Vector & z) {
        for(Natural i=0;i<x.size();i++) 
            z[i]=x[i]*y[i];
    }

    // Identity element, x <- e such that x o e = x.
    static void id(Vector & x) {
        for(Natural i=0;i<x.size();i++) 
            x[i]=Real(1.);
    }
    
    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
    static void linv(Vector const & x,Vector const & y,Vector & z) {
        for(Natural i=0;i<x.size();i++) 
            z[i]=y[i]/x[i];
    }

    // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
    static Real barr(Vector const & x) {
        Real z=Real(0.);
        for(Natural i=0;i<x.size();i++)
            z+=log(x[i]);
        return z;
    }

    // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
    // where y > 0. 
    static Real srch(Vector const & x,Vector const & y) {
        // Line search parameter
        Real alpha=std::numeric_limits <Real>::infinity();

        // Search for the optimal linesearch parameter.
        for(Natural i=0;i<x.size();i++) {
            if(x[i] < Real(0.)) {
                Real alpha0 = -y[i]/x[i];
                alpha = alpha0 < alpha ? alpha0 : alpha;
            }
        }

        return alpha;
    }

    // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
    // operator.
    static void symm(Vector & x) { }
};
//---VectorSpace1---

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosenbrock : public Optizelle::ScalarValuedFunction <double,MyVS> {
    typedef MyVS <double> X;

    // Evaluation of the Rosenbrock function
    double eval(const X::Vector & x) const {
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]));
    }

    // Gradient
    void grad(
        const X::Vector & x,
        X::Vector & g
    ) const {
        g[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0]);
        g[1]=200*(x[1]-sq(x[0]));
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector & x,
        const X::Vector & dx,
        X::Vector & H_dx
    ) const {
    	H_dx[0]= (1200*sq(x[0])-400*x[1]+2)*dx[0]-400*x[0]*dx[1];
        H_dx[1]= -400*x[0]*dx[0] + 200*dx[1];
    }
};

// Define a perfect preconditioner for the Hessian
struct RosenHInv : public Optizelle::Operator <double,MyVS,MyVS> {
public:
    typedef MyVS <double> X;
    typedef X::Vector X_Vector;
private:
    X_Vector& x;
public:
    RosenHInv(X::Vector& x_) : x(x_) {}
    void eval(const X_Vector& dx,X_Vector &result) const {
        double one_over_det=1./(80000.*sq(x[0])-80000.*x[1]+400.);
        result[0]=one_over_det*(200.*dx[0]+400.*x[0]*dx[1]);
        result[1]=one_over_det*
            (400.*x[0]*dx[0]+(1200.*x[0]*x[0]-400.*x[1]+2.)*dx[1]);
    }
};

//---Messaging0---
// Define a custom messaging object
struct MyMessaging : public Optizelle::Messaging {
    // Prints a message
    void print(std::string const & msg) const {
        std::cout << "PRINT:  " << msg << std::endl;
    }

    // Prints an error
    void error(std::string const & msg) const {
        std::cerr << "ERROR:  " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
};
//---Messaging1---

//---Serialization0---
// Define serialization routines for MyVS
namespace Optizelle {
    namespace json {
        template <>
        struct Serialization <double,MyVS> {
            static std::string serialize(
                typename MyVS <double>::Vector const & x
            ) {
                // Create a string with the format 
                // [ x1, x2, ..., xm ].
                std::stringstream x_json;
                x_json.setf(std::ios::scientific);
                x_json.precision(16);
                x_json << "[ ";
                for(Natural i=0;i<x.size()-1;i++)
                    x_json << x[i] << ", ";
                x_json << x.back() << " ]";

                // Return the string
                return x_json.str();
            }
            static MyVS <double>::Vector deserialize(
                typename MyVS <double>::Vector const & x_,
                std::string const & x_json_
            ) {
                // Make a copy of x_json_
                std::string x_json(x_json_);

                // Filter out the commas and brackets from the string
                char formatting[] = "[],";
                for(Natural i=0;i<3;i++) 
                    x_json.erase(
                        std::remove(x_json.begin(),x_json.end(),formatting[i]),
                        x_json.end());

                // Create a new vector that we eventually return
                std::vector <double> x(x_.size());

                // Create a stream out of x_json
                std::stringstream ss(x_json);

                // Read in each of the elements
                for(Natural i=0;i<x.size();i++)
                    ss >> x[i];

                // Return the result 
                return std::move(x);
            }
        };
    }
}
//---Serialization1---

//---RestartManipulator0---
// Define a state manipulator that writes out the optimization state at
// each iteration.
struct MyRestartManipulator
    : Optizelle::StateManipulator <Optizelle::Unconstrained <double,MyVS> >
{
    void eval(
        typename Optizelle::Unconstrained <double,MyVS>
            ::Functions::t const & fns,
        typename Optizelle::Unconstrained <double,MyVS>
            ::State::t & state,
        Optizelle::OptimizationLocation::t const & loc
    ) const {
        switch(loc) {
        // At the end of the optimization iteration, write the restart file
        case Optizelle::OptimizationLocation::EndOfOptimizationIteration: {
            // Create a reasonable file name
            std::stringstream ss;
            ss << "rosenbrock_advanced_api_";
            ss << std::setw(4) << std::setfill('0') << state.iter;
            ss << ".json";

            // Write the restart file
            Optizelle::json::Unconstrained <double,MyVS>::write_restart(
                MyMessaging(),ss.str(),state);
            break;
        } default:
            break;
        }
    }
};
//---RestartManipulator1---

int main(int argc,char* argv[]) {
    // Read in the name for the parameters and optional restart file 
    if(!(argc==2 || argc==3)) {
        std::cerr << "rosenbrock_advanced_api <parameters>" << std::endl;
        std::cerr << "rosenbrock_advanced_api <parameters> <restart>"
            << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string pname(argv[1]);
    std::string rname(argc==3 ? argv[2] : ""); 

    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=-1.2; x[1]=1.;

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <double,MyVS>::State::t state(x);
   
    //---ReadRestart0---
    // If we have a restart file, read in the parameters 
    if(argc==3)
        Optizelle::json::Unconstrained <double,MyVS>::read_restart(
            MyMessaging(),rname,x,state);

    // Read additional parameters from file
    Optizelle::json::Unconstrained <double,MyVS>
        ::read(MyMessaging(),pname,state);
    //---ReadRestart1---

    // Create the bundle of functions 
    Optizelle::Unconstrained <double,MyVS>::Functions::t fns;
    fns.f.reset(new Rosenbrock);
    fns.PH.reset(new RosenHInv(state.x));
    
    //---Solver0---
    // Solve the optimization problem
    Optizelle::Unconstrained <double,MyVS>::Algorithms
        ::getMin(MyMessaging(),fns,state,MyRestartManipulator());
    //---Solver1---

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;

    //---WriteRestart0---
    // Write out the final answer to file
    Optizelle::json::Unconstrained <double,MyVS>::write_restart(
        MyMessaging(),"solution.json",state);
    //---WriteRestart1---
}
