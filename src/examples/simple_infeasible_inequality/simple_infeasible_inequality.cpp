// Optimize a simple optimization problem with an optimal solution of (1/3,1/3),
//
// min x + y
// st  x + 2y >= 1
//     2x + y >= 1
//
// Now, in the case we don't have a starting feasible solution, we can play
// a reformulation trick that adds two scalar variables and allows us to find
// a strictly feasible solution.  Namely,
//
// min x + y
// st  x + 2y >= 1 - z
//     2x + y >= 1 - z
//     epsilon >= w             
//     z = w
//
// Note, most of the time, we're much better off just adding slack variables.
// Basically, this trick is only worthwhile when we don't have a linear system
// solver for the equality constraints added from the slacks since this method
// only adds a single equality constraint.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Define an objective where 
// 
// f(x,y,z,w)=x+y
//
template <typename Real>
struct MyObj : public Optizelle::ScalarValuedFunction <Real,Optizelle::Rm> {
    typedef Optizelle::Rm <Real> X;
    typedef typename X::Vector X_Vector;

    // Evaluation 
    double eval(const X_Vector& x) const {
        return x[0]+x[1];
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad
    ) const {
        grad[0]=Real(1.);
        grad[1]=Real(1.);
        grad[2]=Real(0.);
        grad[3]=Real(0.);
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

// Define a single equality where
//
// g(x,y,z,w) = z - w = 0
//
template <typename Real>
struct MyEq :
    public Optizelle::VectorValuedFunction <Real,Optizelle::Rm,Optizelle::Rm>
{
private:
    Real epsilon;
public:
    typedef Optizelle::Rm <Real> X;
    typedef typename X::Vector X_Vector;
    typedef Optizelle::Rm <Real> Y;
    typedef typename Y::Vector Y_Vector;

    // y=g(x) 
    void eval(
        const X_Vector& x,
        Y_Vector& y
    ) const {
        y[0] = x[2]-x[3];
    }

    // y=g'(x)dx
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Y_Vector& y
    ) const {
        y[0]= dx[2]-dx[3];
    }

    // xhat=g'(x)*dy
    void ps(
        const X_Vector& x,
        const Y_Vector& dy,
        X_Vector& xhat 
    ) const {
        xhat[0]= Real(0.);
        xhat[1]= Real(0.);
        xhat[2]= dy[0];
        xhat[3]= -dy[0];
    }

    // xhat=(g''(x)dx)*dy
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Y_Vector& dy,
        X_Vector& xhat 
    ) const {
        X::zero(xhat);
    }
};


// Define some inequalities where
//
// h(x,y,z,w) = [ x + 2y >= 1 - z  ] 
//              [ 2x + y >= 1 - z  ] 
//              [ epsilon >= w     ]
//
template <typename Real>
struct MyIneq :
    public Optizelle::VectorValuedFunction <Real,Optizelle::Rm,Optizelle::Rm>
{
private:
    Real const & epsilon;
public:
    typedef Optizelle::Rm <Real> X;
    typedef typename X::Vector X_Vector;
    typedef Optizelle::Rm <Real> Z;
    typedef typename Z::Vector Z_Vector;

    // Read in the amount of allowable infeasibility into the problem
    MyIneq(Real const & epsilon_) : epsilon(epsilon_) {}

    // z=h(x) 
    void eval(
        const X_Vector& x,
        Z_Vector& z
    ) const {
        z[0]=x[0]+Real(2.)*x[1]+x[2]-Real(1.);
        z[1]=Real(2.)*x[0]+x[1]+x[2]-Real(1.);
        z[2]=epsilon-x[3];
    }

    // z=h'(x)dx
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Z_Vector& z
    ) const {
        z[0]= dx[0]+Real(2.)*dx[1]+dx[2];
        z[1]= Real(2.)*dx[0]+dx[1]+dx[2];
        z[2]= -dx[3]; 
    }

    // xhat=h'(x)*dz
    void ps(
        const X_Vector& x,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        xhat[0]= dz[0]+Real(2.)*dz[1];
        xhat[1]= Real(2.)*dz[0]+dz[1];
        xhat[2]= dz[0]+dz[1]; 
        xhat[3]= -dz[2];
    }

    // xhat=(h''(x)dx)*dz
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        X::zero(xhat);
    }
};

int main(int argc,char* argv[]){
    // Create a type shortcut
    using Optizelle::Rm;
    typedef double Real;

    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_infeasible_inequality <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Set the amount of infeasibility that we want to allow
    Real const epsilon(1e-8);

    // Generate an initial guess for the primal
    std::vector <Real> x(4);
    x[0]=Real(0.); x[1]=Real(0.); x[2]=Real(5.); x[3]=-Real(5.);

    // Generate a vector for the equality multiplier 
    std::vector <Real> y(1);

    // Generate a vector for the inequality multiplier
    std::vector <Real> z(3);

    // Create an optimization state
    Optizelle::Constrained <Real,Rm,Rm,Rm>::State::t state(x,y,z);

    // Read the parameters from file
    Optizelle::json::Constrained <Real,Rm,Rm,Rm>::read(
        Optizelle::Messaging(),fname,state);
    
    // Create a bundle of functions
    Optizelle::Constrained <Real,Rm,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj <Real>);
    fns.g.reset(new MyEq <Real>);
    fns.h.reset(new MyIneq <Real>(epsilon));

    // Solve the optimization problem
    std::cout << std::endl << "Solving the optimization problem." << std::endl;
    Optizelle::Constrained <Real,Rm,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << std::scientific << std::setprecision(16)
        << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::Constrained <Real,Rm,Rm,Rm>::write_restart(
        Optizelle::Messaging(),"solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
