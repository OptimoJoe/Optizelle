// Optimize a simple optimization problem with an optimal 
// solution of (1/3,1/3)

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

// Grab Optizelle's Natural type
using Optizelle::Natural;

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
    static void scal(Real const & alpha, Vector & x) {
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
    static void axpy(Real const & alpha, Vector const & x, Vector & y) {
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

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Define a simple objective where 
// 
// f(x,y)=(x+1)^2+(y+1)^2
//
struct MyObj : public Optizelle::ScalarValuedFunction <double,MyVS> {
    typedef MyVS <double> X;

    // Evaluation 
    double eval(const X::Vector& x) const {
        return sq(x[0]+1.)+sq(x[1]+1.);
    }

    // Gradient
    void grad(
        X::Vector const & x,
        X::Vector & grad
    ) const {
        grad[0]=2*x[0]+2;
        grad[1]=2*x[1]+2;
    }

    // Hessian-vector product
    void hessvec(
        X::Vector const & x,
        X::Vector const & dx,
        X::Vector & H_dx
    ) const {
        H_dx[0]=2.*dx[0]; 
        H_dx[1]=2.*dx[1]; 
    }
};

// Define a simple equality
//
// g(x,y)= [ x + 2y = 1 ] 
//
struct MyEq :public Optizelle::VectorValuedFunction<double,MyVS,MyVS> {
    typedef MyVS <double> X;
    typedef MyVS <double> Y;

    // y=g(x) 
    void eval(
        X::Vector const & x,
        Y::Vector & y
    ) const {
        y[0]=x[0]+2.*x[1]-1.;
    }

    // y=g'(x)dx
    void p(
        X::Vector const & x,
        X::Vector const & dx,
        Y::Vector & y
    ) const {
        y[0]= dx[0]+2.*dx[1];
    }

    // xhat=g'(x)*dy
    void ps(
        X::Vector const & x,
        Y::Vector const & dy,
        X::Vector & xhat 
    ) const {
        xhat[0]= dy[0];
        xhat[1]= 2.*dy[0];
    }

    // xhat=(g''(x)dx)*dy
    void pps(
        X::Vector const & x,
        X::Vector const & dx,
        Y::Vector const & dy,
        X::Vector & xhat 
    ) const {
        X::zero(xhat);
    }
};

// Define a simple inequality
//
// h(x,y)= [ 2x + y >= 1 ] 
//
struct MyIneq :public Optizelle::VectorValuedFunction<double,MyVS,MyVS> {
    typedef MyVS <double> X;
    typedef MyVS <double> Z;

    // z=h(x) 
    void eval(
        X::Vector const & x,
        Z::Vector & z
    ) const {
        z[0]=2.*x[0]+x[1]-1.;
    }

    // z=h'(x)dx
    void p(
        X::Vector const & x,
        X::Vector const & dx,
        Z::Vector & z
    ) const {
        z[0]= 2.*dx[0]+dx[1];
    }

    // xhat=h'(x)*dz
    void ps(
        X::Vector const & x,
        Z::Vector const & dz,
        X::Vector & xhat 
    ) const {
        xhat[0]= 2.*dz[0];
        xhat[1]= dz[0];
    }

    // xhat=(h''(x)dx)*dz
    void pps(
        X::Vector const & x,
        X::Vector const & dx,
        Z::Vector const & dz,
        X::Vector & xhat 
    ) const {
        X::zero(xhat);
    }
};

//---Serialization0---
// Define serialization routines for MyVS
namespace Optizelle {
    namespace json {
        template <>
        struct Serialization <double,MyVS> {
            static std::string serialize(
                typename MyVS <double>::Vector const & x,
                std::string const & name,
                Natural const & iter
            ) {
                // Create the filename where we put our vector
                std::stringstream fname;
                fname << "./restart/";
                fname << name << ".";
                fname << std::setw(4) << std::setfill('0') << iter;
                fname << ".txt";

                // Actually write the vector there
                std::ofstream fout(fname.str());
                if(fout.fail()) {
                    std::stringstream msg;
                    msg << "While writing the variable " << name 
                        << " to file on iteration " << iter
                        << ", unable to open the file: "
                        << fname.str() << ".";
                    throw Optizelle::Exception::t(msg.str());
                }
                fout.setf(std::ios::scientific);
                fout.precision(16);
                for(Natural i=0;i<x.size();i++)
                    fout << x[i] << std::endl;

                // Close out the file
                fout.close();

                // Use this filename as the json string 
                std::stringstream x_json;
                x_json << "\"" << fname.str() << "\"";
                return x_json.str();
            }
            static MyVS <double>::Vector deserialize(
                typename MyVS <double>::Vector const & x_,
                std::string const & x_json_
            ) {
                // Make a copy of x_json_
                auto x_json = x_json_;

                // Filter out the quotes and newlines from the string
                auto formatting = "\"\n";
                for(auto i=0;i<2;i++) 
                    x_json.erase(
                        std::remove(x_json.begin(),x_json.end(),formatting[i]),
                        x_json.end());

                // Open the file for reading 
                std::ifstream fin(x_json.c_str());
                if(!fin.is_open())
                    throw Optizelle::Exception::t(
                        "Error while opening the file " + x_json + ": " +
                        strerror(errno));

                // Create a new vector that we eventually return
                auto x = std::vector <double> (x_.size());

                // Read in each of the elements
                for(auto i=0;i<x.size();i++)
                    fin >> x[i];

                // Return the result 
                return std::move(x);
            }
        };
    }
}
//---Serialization1---

// Define a state manipulator that writes out the optimization state at
// each iteration.
struct MyRestartManipulator : Optizelle::StateManipulator <
    Optizelle::Constrained <double,MyVS,MyVS,MyVS> >
{
    void eval(
        typename Optizelle::Constrained <double,MyVS,MyVS,MyVS>
            ::Functions::t const & fns,
        typename Optizelle::Constrained <double,MyVS,MyVS,MyVS>
            ::State::t & state,
        Optizelle::OptimizationLocation::t const & loc
    ) const {
        switch(loc) {
        // At the end of the optimization iteration, write the restart file
        case Optizelle::OptimizationLocation::EndOfOptimizationIteration: {
            // Create a reasonable file name
            std::stringstream ss;
            ss << "simple_constrained_advanced_api_";
            ss << std::setw(4) << std::setfill('0') << state.iter;
            ss << ".json";

            // Write the restart file
            Optizelle::json::Constrained <double,MyVS,MyVS,MyVS>::write_restart(
                ss.str(),state);
            break;
        } default:
            break;
        }
    }
};

int main(int argc,char* argv[]){
    // Read in the name for the parameters and optional restart file 
    if(!(argc==2 || argc==3)) {
        std::cerr << "simple_constrained_advanced_api <parameters>"<< std::endl;
        std::cerr << "simple_constrained_advanced_api <parameters> <restart>"
            << std::endl;
        exit(EXIT_FAILURE);
    }
    auto pname = argv[1];
    auto rname = argc==3 ? argv[2] : ""; 

    // Generate an initial guess for the primal
    auto x = std::vector <double> {2.1, 1.1};

    // Allocate memory for equality multiplier 
    auto y = std::vector <double> (1);

    // Allocate memory for the inequality multiplier 
    auto z = std::vector <double> (1);

    // Create an optimization state
    Optizelle::Constrained <double,MyVS,MyVS,MyVS>::State::t state(x,y,z);
    
    // If we have a restart file, read in the parameters 
    if(argc==3)
        Optizelle::json::Constrained <double,MyVS,MyVS,MyVS>::read_restart(
            rname,x,y,z,state);

    // Read the parameters from file
    Optizelle::json::Constrained <double,MyVS,MyVS,MyVS>::read(pname,state);
    
    // Create a bundle of functions
    Optizelle::Constrained <double,MyVS,MyVS,MyVS>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.g.reset(new MyEq);
    fns.h.reset(new MyIneq);
    
    // Solve the optimization problem
    Optizelle::Constrained <double,MyVS,MyVS,MyVS>::Algorithms
        ::getMin(Optizelle::Messaging::stdout,fns,state,MyRestartManipulator());

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    std::cout << std::scientific << std::setprecision(16)
        << "The optimal point is: (" << state.x[0] << ','
        << state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::Constrained <double,MyVS,MyVS,MyVS>::write_restart(
        "solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
