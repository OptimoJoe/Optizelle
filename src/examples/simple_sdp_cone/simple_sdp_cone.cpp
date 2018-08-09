// Optimize a simple problem with an optimal solution of (0.5,.25)

#include <iostream>
#include <iomanip>
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"

// Create some type shortcuts
using Optizelle::Rm;
using Optizelle::SQL;
typedef double Real;

// Define a simple objective where
//
// f(x,y)=-x+y
//
struct MyObj : public Optizelle::ScalarValuedFunction <Real,Rm> {
    typedef Rm <Real> X;

    // Evaluation
    double eval(X::Vector const & x) const {
        return -x[0]+x[1];
    }

    // Gradient
    void grad(
        X::Vector const & x,
        X::Vector & grad
    ) const {
        grad[0]=Real(-1.);
        grad[1]=Real(1.);
    }

    // Hessian-vector product
    void hessvec(
        X::Vector const & x,
        X::Vector const & dx,
        X::Vector & H_dx
    ) const {
        H_dx[0]= Real(0.);
        H_dx[1]= Real(0.);
    }
};

// Define a simple SDP inequality
//
// h(x,y) = [ y x ] >= 0
//          [ x 1 ]
//
struct MyIneq : public Optizelle::VectorValuedFunction <Real,Rm,SQL> {
    typedef Rm <Real> X;
    typedef SQL <Real> Z;

    // z=h(x)
    void eval(
        X::Vector const & x,
        Z::Vector & z
    ) const {
        z(1,1,1)=x[1];
        z(1,1,2)=x[0];
        z(1,2,1)=x[0];
        z(1,2,2)=Real(1.);
    }

    // z=h'(x)dx
    void p(
        X::Vector const & x,
        X::Vector const & dx,
        Z::Vector & z
    ) const {
        z(1,1,1)=dx[1];
        z(1,1,2)=dx[0];
        z(1,2,1)=dx[0];
        z(1,2,2)=Real(0.);
    }

    // xhat=h'(x)*dz
    void ps(
        X::Vector const & x,
        Z::Vector const & dz,
        X::Vector & xhat
    ) const {
        xhat[0]= dz(1,1,2)+dz(1,2,1);
        xhat[1]= dz(1,1,1);
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

int main(int argc,char* argv[]){
    // Create some type shortcuts
    typedef Rm <Real>::Vector Rm_Vector;
    typedef SQL <Real>::Vector SQL_Vector;

    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_sdp_cone <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    auto fname = argv[1];

    // Generate an initial guess for the primal
    auto x = Rm_Vector({1.2,3.1});

    // Allocate memory for the dual
    auto z = SQL_Vector ({Optizelle::Cone::Semidefinite},{2});

    // Create an optimization state
    Optizelle::InequalityConstrained <Real,Rm,SQL>::State::t state(x,z);

    // Read the parameters from file
    Optizelle::json::InequalityConstrained <Real,Rm,SQL>::read(fname,state);

    // Create a bundle of functions
    Optizelle::InequalityConstrained <Real,Rm,SQL>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Solve the optimization problem
    Optizelle::InequalityConstrained <Real,Rm,SQL>
        ::Algorithms::getMin(Optizelle::Messaging::stdout,fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << std::setprecision(16) << std::scientific
        << "The optimal point is: (" << state.x[0] << ','
        << state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::InequalityConstrained <Real,Rm,SQL>
        ::write_restart("solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
