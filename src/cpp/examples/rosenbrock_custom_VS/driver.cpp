// This example minimizes the Rosenbrock function from scratch.  Meaning,
// it runs through a complete example from defining a valid vector space,
// to setting parameters, to solving the problem.  For reference, the optimal
// solution to the Rosenbrock function is (1,1).

#include <vector>
#include <iostream>
#include <string>
#include "optizelle/optizelle.h"

// Grab Optizelle's Natural type
using Optizelle::Natural;

// Defines the vector space used for optimization.
template <typename Real>
struct MyHS { 
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
};

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosen : public Optizelle::ScalarValuedFunction <double,MyHS> {
    typedef MyHS <double> X;

    // Evaluation of the Rosenbrock function
    double operator () (const X::Vector & x) const {
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

int main(){

    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=-1.2; x[1]=1.;

    // Create a direction for the finite difference tests
    std::vector <double> dx(2);
    dx[0]=-.5; dx[1]=.5;
    
    // Create another direction for the finite difference tests
    std::vector <double> dxx(2);
    dxx[0]=.75; dxx[1]=.25;

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <double,MyHS>::State::t state(x);

    // Setup some algorithmic parameters
    #if 1
    // Trust-Region Newton's method
    state.H_type = Optizelle::Operators::UserDefined;
    state.iter_max = 50;
    state.eps_krylov = 1e-10;
    #endif

    // BFGS
    #if 0
    state.algorithm_class = Optizelle::AlgorithmClass::LineSearch;
    state.dir = Optizelle::LineSearchDirection::BFGS;
    state.stored_history = 10;
    state.iter_max = 100;
    #endif
    
    // Newton-CG 
    #if 0
    state.algorithm_class = Optizelle::AlgorithmClass::LineSearch;
    state.dir = Optizelle::LineSearchDirection::NewtonCG;
    state.H_type = Optizelle::Operators::UserDefined;
    state.eps_krylov = 1e-2;
    state.iter_max = 50;
    #endif

    // Create the bundle of functions 
    Optizelle::Unconstrained <double,MyHS>::Functions::t fns;
    fns.f.reset(new Rosen);
    
    // Solve the optimization problem
    Optizelle::Unconstrained <double,MyHS>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;
}
