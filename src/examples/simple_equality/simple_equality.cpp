// Optimize a simple optimization problem with an optimal solution 
// of (2-sqrt(2)/2,2-sqrt(2)/2).

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Define a simple objective where 
// 
// f(x,y)=x^2+y^2
//
struct MyObj
    : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;

    // Evaluation 
    double eval(const X::Vector& x) const {
        return sq(x[0])+sq(x[1]);
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& grad
    ) const {
        grad[0]=2.*x[0];
        grad[1]=2.*x[1];
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
        H_dx[0]=2.*dx[0]; 
        H_dx[1]=2.*dx[1]; 
    }
};

//---EqualityConstraint0---
// Define a simple equality constraint
//
// g(x,y)= [ (x-2)^2 + (y-2)^2 = 1 ] 
//
struct MyEq
    :public Optizelle::VectorValuedFunction<double,Optizelle::Rm,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;

    // y=g(x) 
    void eval(
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[0] = sq(x[0]-2.)+sq(x[1]-2.)-1.;
    }

    // y=g'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y[0] = 2.*(x[0]-2.)*dx[0]+2.*(x[1]-2.)*dx[1];
    }

    // z=g'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0] = 2.*(x[0]-2.)*dy[0];
        z[1] = 2.*(x[1]-2.)*dy[0];
    }

    // z=(g''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0] = 2.*dx[0]*dy[0];
        z[1] = 2.*dx[1]*dy[0];
    }
};
//---EqualityConstraint1---

int main(int argc,char* argv[]){
    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_equalty <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Create a type shortcut
    using Optizelle::Rm;

    // Generate an initial guess 
    std::vector <double> x(2);
    x[0]=2.1; x[1]=1.1;

    // Allocate memory for the equality multiplier 
    std::vector <double> y(1);

    // Create an optimization state
    Optizelle::EqualityConstrained <double,Rm,Rm>::State::t state(x,y);

    // Read the parameters from file
    Optizelle::json::EqualityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::read(Optizelle::Messaging(),fname,state);
    
    // Create a bundle of functions
    Optizelle::EqualityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.g.reset(new MyEq);

    // Solve the optimization problem
    Optizelle::EqualityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    std::cout << std::scientific << std::setprecision(16)
        << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::EqualityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::write_restart(Optizelle::Messaging(),"solution.json",state);

    // Return that we've exited successfuly
    return EXIT_SUCCESS;
}
