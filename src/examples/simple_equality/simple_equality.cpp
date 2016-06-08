// Optimize a simple optimization problem with an optimal solution 
// of (2-sqrt(2)/2,2-sqrt(2)/2).

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

//---Objective0---
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
    double eval(X::Vector const & x) const {
        return sq(x[0])+sq(x[1]);
    }

    // Gradient
    void grad(
        X::Vector const & x,
        X::Vector & grad
    ) const {
        grad[0]=2.*x[0];
        grad[1]=2.*x[1];
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
//---Objective1---

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
        X::Vector const & x,
        Y::Vector & y
    ) const {
        y[0] = sq(x[0]-2.)+sq(x[1]-2.)-1.;
    }

    // y=g'(x)dx
    void p(
        X::Vector const & x,
        X::Vector const & dx,
        Y::Vector & y
    ) const {
        y[0] = 2.*(x[0]-2.)*dx[0]+2.*(x[1]-2.)*dx[1];
    }

    // xhat=g'(x)*dy
    void ps(
        X::Vector const & x,
        Y::Vector const & dy,
        X::Vector & xhat 
    ) const {
        xhat[0] = 2.*(x[0]-2.)*dy[0];
        xhat[1] = 2.*(x[1]-2.)*dy[0];
    }

    // xhat=(g''(x)dx)*dy
    void pps(
        X::Vector const & x,
        X::Vector const & dx,
        Y::Vector const & dy,
        X::Vector & xhat 
    ) const {
        xhat[0] = 2.*dx[0]*dy[0];
        xhat[1] = 2.*dx[1]*dy[0];
    }
};
//---EqualityConstraint1---

//---Preconditioner0---
// Define a Schur preconditioner for the equality constraints 
struct MyPrecon:
    public Optizelle::Operator <double,Optizelle::Rm,Optizelle::Rm>
{
public:
    typedef Optizelle::Rm <double> X;
    typedef X::Vector X_Vector;
    typedef Optizelle::Rm <double> Y;
    typedef Y::Vector Y_Vector;
private:
    X_Vector& x;
public:
    MyPrecon(X::Vector& x_) : x(x_) {}
    void eval(Y_Vector const & dy,Y_Vector & result) const {
        result[0]=dy[0]/sq(4.*(x[0]-2.)+4.*sq(x[1]-2.));
    }
};
//---Preconditioner1---

int main(int argc,char* argv[]){
    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_equalty <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    auto fname = argv[1];

    // Create a type shortcut
    using Optizelle::Rm;

    //---State0---
    // Generate an initial guess 
    auto x = std::vector <double> {2.1, 1.1};

    // Allocate memory for the equality multiplier 
    auto y = std::vector <double> (1);

    // Create an optimization state
    Optizelle::EqualityConstrained <double,Rm,Rm>::State::t state(x,y);
    //---State1---

    //---Parameters0---
    // Read the parameters from file
    Optizelle::json::EqualityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::read(fname,state);
    //---Parameters1---
   
    //---Functions0---
    // Create a bundle of functions
    Optizelle::EqualityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.g.reset(new MyEq);
    fns.PSchur_left.reset(new MyPrecon(state.x));
    //---Functions1---

    //---Solver0---
    // Solve the optimization problem
    Optizelle::EqualityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging::stdout,fns,state);
    //---Solver1---

    //---Extract0---
    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    std::cout << std::scientific << std::setprecision(16)
        << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;
    //---Extract1---

    // Write out the final answer to file
    Optizelle::json::EqualityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::write_restart("solution.json",state);

    // Return that we've exited successfuly
    return EXIT_SUCCESS;
}
