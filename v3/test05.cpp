#include <vector>
#include <iostream>
#include <string>
#include "peopt.h"
#include "vspaces.h"

// Define a simple optimization problem with an optimal solution of (2.5,2.5)

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define a simple objective where 
// 
// f(x,y)=(x-3)^2+(y-2)^2
//
struct MyObj : public peopt::ScalarValuedFunction <double,peopt::Rm> {
    typedef double Real;
    typedef peopt::Rm <Real> X;

    // Evaluation 
    double operator () (const X::Vector& x) const {
        return sq(x[0]-Real(3.))+sq(x[1]-Real(2.));
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& g
    ) const {
        g[0]=2*x[0]-6;
        g[1]=2*x[1]-4;
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
    	H_dx[0]= Real(2.)*dx[0]; 
        H_dx[1]= Real(2.)*dx[1]; 
    }
};

// Define a simple SOCP inequality 
//
// g(x,y) = [ y >= |x| ] 
// g(x,y) =  (y,x) >=_Q 0
//
struct MyIneq:public peopt::VectorValuedFunction <double,peopt::Rm,peopt::SQL> {
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Y;
    typedef double Real;

    // y=f(x) 
    void operator () (
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y(1,1)=x[1];
        y(1,2)=x[0];
    }

    // y=f'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y(1,1)= dx[1];
        y(1,2)= dx[0];
    }

    // z=f'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0]= dy(1,2);
        z[1]= dy(1,1);
    }

    // z=(f''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        X::zero(z);
    }
};

int main(){
    // Create some type shortcuts
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <unsigned int> sizes(1); sizes[0]=2;
    std::vector <peopt::Cone::t> types(1); types[0]=peopt::Cone::Quadratic;
    Z_Vector z(peopt::Messaging(),types,sizes);
    Z::id(z);

    // Create a direction for the finite difference tests
    X_Vector dx(2);
    dx[0]=-.5; dx[1]=.5;
    
    // Create another direction for the finite difference tests
    X_Vector dxx(2);
    dxx[0]=.75; dxx[1]=.25;
    
    // Create a vector in the codomain of the function used in the inequality
    // constraints 
    Z_Vector dz;
    Z::init(z,dz);
    dz(1,1)=.3; dz(1,2)=-.12; 
    
    // Create an optimization state
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);
    
    // Do some finite difference tests on the objective
    peopt::Diagnostics::gradientCheck <>
        (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianCheck <>
        (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianSymmetryCheck <>
        (peopt::Messaging(),*(fns.f),x,dx,dxx);
    
    // Do some finite difference tests on the constraint 
    peopt::Diagnostics::derivativeCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dz);
    peopt::Diagnostics::derivativeAdjointCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dz);
    peopt::Diagnostics::secondDerivativeCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dz);

    // Setup the optimization problem
    #if 1
    // Newton's method
    state.H_type = peopt::Operators::External;
    state.iter_max = 100;
    state.eps_krylov = 1e-10;
    state.eps_s = 1e-16;
    state.eps_g = 1e-10;
    state.sigma = 0.10;
    state.gamma = 0.95;
    #endif

    // BFGS
    #if 0
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::BFGS;
    state.stored_history = 10;
    state.iter_max = 300;
    state.sigma = 0.10;
    state.gamma = 0.95;
    state.eps_s = 1e-16;
    #endif
    
    // Newton-CG 
    #if 0
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::NewtonCG;
    state.H_type = peopt::Operators::External;
    state.eps_krylov = 1e-10;
    state.iter_max = 300;
    state.eps_s = 1e-16;
    state.eps_g = 1e-10;
    state.sigma = 0.10;
    state.gamma = 0.95;
    #endif

    // Solve the optimization problem
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        peopt::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    const X_Vector& opt_x=*(state.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;
}
