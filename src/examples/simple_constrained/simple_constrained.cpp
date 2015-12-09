// Optimize a simple optimization problem with an optimal 
// solution of (1/3,1/3)

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
// f(x,y)=(x+1)^2+(y+1)^2
//
struct MyObj
    : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;

    // Evaluation 
    double eval(const X::Vector& x) const {
        return sq(x[0]+1.)+sq(x[1]+1.);
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& g
    ) const {
        g[0]=2*x[0]+2;
        g[1]=2*x[1]+2;
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

// Define a simple equality
//
// g(x,y)= [ x + 2y = 1 ] 
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
        y[0]=x[0]+2.*x[1]-1.;
    }

    // y=g'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y[0]= dx[0]+2.*dx[1];
    }

    // z=g'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0]= dy[0];
        z[1]= 2.*dy[0];
    }

    // z=(g''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        X::zero(z);
    }
};

// Define a simple inequality
//
// h(x,y)= [ 2x + y >= 1 ] 
//
struct MyIneq
    :public Optizelle::VectorValuedFunction<double,Optizelle::Rm,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;

    // y=h(x) 
    void eval(
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[0]=2.*x[0]+x[1]-1.;
    }

    // y=h'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y[0]= 2.*dx[0]+dx[1];
    }

    // z=h'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0]= 2.*dy[0];
        z[1]= dy[0];
    }

    // z=(h''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        X::zero(z);
    }
};

int main(int argc,char* argv[]){
    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_constrained <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Create a type shortcut
    using Optizelle::Rm;

    // Generate an initial guess for the primal
    std::vector <double> x(2);
    x[0]=2.1; x[1]=1.1;

    // Allocate memory for equality multiplier 
    std::vector <double> y(1);

    // Allocate memory for the inequality multiplier 
    std::vector <double> z(1);

    // Create an optimization state
    Optizelle::Constrained <double,Rm,Rm,Rm>::State::t state(x,y,z);

    // Read the parameters from file
    Optizelle::json::Constrained <double,Rm,Rm,Rm>::read(
        Optizelle::Messaging(),fname,state);
    
    // Create a bundle of functions
    Optizelle::Constrained <double,Rm,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.g.reset(new MyEq);
    fns.h.reset(new MyIneq);
    
    // Solve the optimization problem
    Optizelle::Constrained <double,Rm,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    std::cout << std::scientific << std::setprecision(16)
        << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::Constrained <double,Rm,Rm,Rm>::write_restart(
        Optizelle::Messaging(),"solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
