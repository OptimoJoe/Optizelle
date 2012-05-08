#include <vector>
#include <iostream>
#include <string>
#include "peopt.h"

// Defines the vector space used for optimization.
template <typename Real>
struct MyHS { 
    typedef std::vector <Real> Vector;

    // Create an empty, uninitialized vector
    static Vector create() {
        return std::vector <Real> ();
    }
    
    // Memory allocation and size setting
    static void init(const Vector& x, Vector& y) {
        y.resize(x.size());
    }

    // y <- x (Shallow.  No memory allocation.)
    static void copy(const Vector& x, Vector& y) {
        for(unsigned int i=0;i<x.size();i++){
            y[i]=x[i];
        }
    }

    // x <- alpha * x
    static void scal(const Real& alpha, Vector& x) {
        for(unsigned int i=0;i<x.size();i++){
            x[i]=alpha*x[i];
        }
    }

    // x <- 0 
    static void zero(Vector& x) {
        for(unsigned int i=0;i<x.size();i++){
            x[i]=0.;
        }
    }

    // y <- alpha * x + y
    static void axpy(const Real& alpha, const Vector& x, Vector& y) {
        for(unsigned int i=0;i<x.size();i++){
            y[i]=alpha*x[i]+y[i];
        }
    }

    // innr <- <x,y>
    static Real innr(const Vector& x,const Vector& y) {
        Real z=0;
        for(unsigned int i=0;i<x.size();i++)
            z+=x[i]*y[i];
        return z;
    }

    // Jordan product, z <- x o y
    static void prod(const Vector& x, const Vector& y, Vector& z) {
        for(unsigned int i=0;i<x.size();i++) 
            z[i]=x[i]*y[i];
    }

    // Identity element, x <- e such that x o e = x
    static void id(Vector& x) {
        for(unsigned int i=0;i<x.size();i++) 
            x[i]=1.;
    }
    
    // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
    static void linv(const Vector& x,const Vector& y,Vector& z) {
        for(unsigned int i=0;i<x.size();i++) 
            z[i]=(1/x[i])*y[i];
    }

    // Barrier function, barr <- barr(x) where x o grad barr(x) = e
    static Real barr(const Vector& x) {
        Real z=0;
        for(unsigned int i=0;i<x.size();i++)
            z+=std::log(x[i]);
        return z;
    }

    // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
    // where y > 0.  If the argmax is infinity, then return Real(-1.).
    static Real srch(const Vector& x,const Vector& y) {
        Real alpha=Real(-1.);
        for(unsigned int i=0;i<x.size();i++) {
            if(x[i] < 0) {
                Real alpha0;
                alpha0 = -y[i]/x[i];
                if(alpha==Real(-1.) || alpha0 < alpha)
                    alpha=alpha0;
            }
        }
        return alpha;
    }
};

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define a simple objective where 
// 
// f(x,y)=(x+1)^2+(y+1)^2
//
struct MyObj : public peopt::ScalarValuedFunction <double,MyHS> {
    typedef MyHS <double> X;
    typedef double Real;

    // Evaluation 
    double operator () (const X::Vector& x) const {
        return sq(x[0]+Real(1.))+sq(x[1]+Real(1.));
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& g
    ) const {
        g[0]=2*x[0]+1;
        g[1]=2*x[1]+1;
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

// Define simple inequalities 
//
// g(x,y)= [ x + 2y >= 1 ] 
//         [ 2x + y >= 1 ] 
//
struct MyIneq : public peopt::VectorValuedFunction <double,MyHS,MyHS> {
    typedef MyHS <double> X;
    typedef MyHS <double> Y;
    typedef double Real;

    // y=f(x) 
    void operator () (
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[0]=x[0]+Real(2.)*x[1]-Real(1.);
        y[1]=Real(2.)*x[0]+x[1]-Real(1.);
    }

    // y=f'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y[0]= dx[0]+Real(2.)*dx[1];
        y[1]= Real(2.)*dx[0]+dx[1];
    }

    // z=f'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0]= dy[0]+Real(2.)*dy[1];
        z[1]= Real(2.)*dy[0]+dy[1];
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

    // Generate an initial guess for the primal
    std::vector <double> x(2);
    x[0]=1.1; x[1]=1.1;

    // Generate an initial guess for the dual
    std::vector <double> z(2);
    z[0]=1.; z[1]=1.;

    // Create a direction for the finite difference tests
    std::vector <double> dx(2);
    dx[0]=-.5; dx[1]=.5;
    
    // Create another direction for the finite difference tests
    std::vector <double> dxx(2);
    dxx[0]=.75; dxx[1]=.25;
    
    // Create a vector in the codomain of the function used in the inequality
    // constraints 
    std::vector <double> dy(2);
    dy[0]=.3; dy[1]=-.12; 
    
    // Create an optimization state
    peopt::InequalityConstrained <double,MyHS,MyHS>::State::t state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,MyHS,MyHS>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);
    #if 0
    peopt::InequalityConstrained <double,MyHS,MyHS>::Functions::init(
        peopt::Messaging(),state,fns);
    #endif
    
    // Do some finite difference tests on the objective
    peopt::Diagnostics::gradientCheck <>
        (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianCheck <>
        (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianSymmetryCheck <>
        (peopt::Messaging(),*(fns.f),x,dx,dxx);
    
    // Do some finite difference tests on the constraint 
    peopt::Diagnostics::derivativeCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dy);
    peopt::Diagnostics::derivativeAdjointCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dy);
    peopt::Diagnostics::secondDerivativeCheck <>
        (peopt::Messaging(),*(fns.h),x,dx,dy);

    // Setup the optimization problem
    #if 0
    // Newton's method
    state.H_type = peopt::Operators::External;
    state.iter_max = 100;
    state.eps_krylov = 1e-10;
    state.eps_s = 1e-16;
    state.eps_g = 1e-10;
    state.sigma = 0.25;
    state.gamma = 0.95;
    #endif

    // BFGS
    #if 0
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::BFGS;
    state.stored_history = 10;
    state.iter_max = 300;
    #endif
    
    // Newton-CG 
    #if 1
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::NewtonCG;
    state.H_type = peopt::Operators::External;
    state.eps_krylov = 1e-10;
    state.iter_max = 300;
    state.eps_s = 1e-16;
    state.eps_g = 1e-10;
    state.sigma = 0.75;
    state.gamma = 0.95;
    #endif

    // Solve the optimization problem
    peopt::InequalityConstrained <double,MyHS,MyHS>::Algorithms
        ::getMin(peopt::Messaging(1),fns,state);

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;
}
