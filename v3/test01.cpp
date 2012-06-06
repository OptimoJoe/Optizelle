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

    // norm <- ||x||
    static Real norm(const Vector& x) {
        return sqrt(innr(x,x));
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

// Cubes its input
template <typename Real>
Real cub(Real x){
    return x*x*x; 
}

// Quads its input
template <typename Real>
Real quad(Real x){
    return x*x*x*x; 
}

// Quints its input
template <typename Real>
Real quint(Real x){
    return x*x*x*x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosen : public peopt::ScalarValuedFunction <double,MyHS> {
    typedef MyHS <double> X;

    // Evaluation of the Rosenbrock function
    double operator () (const X::Vector& x) const {
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]));
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& g
    ) const {
        g[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0]);
        g[1]=200*(x[1]-sq(x[0]));
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
    	H_dx[0]= (1200*sq(x[0])-400*x[1]+2)*dx[0]-400*x[0]*dx[1];
        H_dx[1]= -400*x[0]*dx[0] + 200*dx[1];
    }
};

// Define some utility function where
//
// g(x,y)= [ cos(x)sin(y) ]
//         [ 3 x^2 y + y^3]
//         [ log(x) + 3y^5]
//
struct Utility  : public peopt::VectorValuedFunction <double,MyHS,MyHS> {
    typedef MyHS <double> X;
    typedef MyHS <double> Y;

    // y=f(x) 
    void operator () (
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[0]=cos(x[0])*sin(x[1]);
        y[1]=3.*sq(x[0])*x[1]+cub(x[1]);
        y[2]=log(x[0])+3.*quint(x[1]);
    }

    // y=f'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        y[0]= -sin(x[0])*sin(x[1])*dx[0]
              +cos(x[0])*cos(x[1])*dx[1];
        y[1]= 6.*x[0]*x[1]*dx[0]
              +(3.*sq(x[0])+3.*sq(x[1]))*dx[1];
        y[2]= 1./x[0]*dx[0]
              +15.*quad(x[1])*dx[1];
    }

    // z=f'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0]= -sin(x[0])*sin(x[1])*dy[0]
              +6.*x[0]*x[1]*dy[1]
              +1./x[0]*dy[2];
        z[1]= cos(x[0])*cos(x[1])*dy[0]
              +(3.*sq(x[0])+3.*sq(x[1]))*dy[1]
              +15.*quad(x[1])*dy[2];
    }

    // z=(f''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        z[0] = (-cos(x[0])*dx[0]*sin(x[1])-sin(x[0])*cos(x[1])*dx[1])*dy[0]
               +(6.*dx[0]*x[1] + 6.*x[0]*dx[1])*dy[1]
               +(-1./sq(x[0])*dx[0])*dy[2];
        z[1] = (-sin(x[0])*dx[0]*cos(x[1])-cos(x[0])*sin(x[1])*dx[1])*dy[0]
               +(6.*x[0]*dx[0]+6.*x[1]*dx[1])*dy[1]
               +(60.*cub(x[1])*dx[1])*dy[2];
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

    // Create a vector in the codomain of the utility function
    std::vector <double> dy(3);
    dy[0]=.3; dy[1]=-.12; dy[2]=1.2;

    // Create an optimization state
    peopt::Unconstrained <double,MyHS>::State::t ustate(x);
    peopt::EqualityConstrained <double,MyHS,MyHS>::State::t estate(x,x);
    peopt::InequalityConstrained <double,MyHS,MyHS>::State::t  istate(x,x);
    peopt::Constrained <double,MyHS,MyHS,MyHS>::State::t cstate(x,x,x);

    // Do a capture and release of the state
    peopt::Constrained <double,MyHS,MyHS,MyHS>::X_Vectors xs;
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Y_Vectors ys;
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Z_Vectors zs;
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Reals reals;
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Nats nats;
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Params params;

    peopt::Unconstrained <double,MyHS>::Restart
        ::release(ustate,xs,reals,nats,params);
    peopt::Unconstrained <double,MyHS>::Restart
        ::capture(peopt::Messaging(),ustate,xs,reals,nats,params);
   
    xs.first.clear(); ys.first.clear(); zs.first.clear();
    reals.first.clear(); nats.first.clear(); params.first.clear();
    xs.second.clear(); ys.second.clear(); zs.second.clear();
    reals.second.clear(); nats.second.clear(); params.second.clear();
    peopt::EqualityConstrained <double,MyHS,MyHS>::Restart
        ::release(estate,xs,ys,reals,nats,params);
    peopt::EqualityConstrained <double,MyHS,MyHS>::Restart
        ::capture(peopt::Messaging(),estate,xs,ys,reals,nats,params);
    
    xs.first.clear(); ys.first.clear(); zs.first.clear();
    reals.first.clear(); nats.first.clear(); params.first.clear();
    xs.second.clear(); ys.second.clear(); zs.second.clear();
    reals.second.clear(); nats.second.clear(); params.second.clear();
    istate.mu_trg=.25;
    peopt::InequalityConstrained <double,MyHS,MyHS>::Restart
        ::release(istate,xs,zs,reals,nats,params);
    peopt::InequalityConstrained <double,MyHS,MyHS>::Restart
        ::capture(peopt::Messaging(),istate,xs,zs,reals,nats,params);
    
    xs.first.clear(); ys.first.clear(); zs.first.clear();
    reals.first.clear(); nats.first.clear(); params.first.clear();
    xs.second.clear(); ys.second.clear(); zs.second.clear();
    reals.second.clear(); nats.second.clear(); params.second.clear();
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Restart
        ::release(cstate,xs,ys,zs,reals,nats,params);
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Restart
        ::capture(peopt::Messaging(),cstate,xs,ys,zs,reals,nats,params);

    // Create a bundle of functions
    peopt::Unconstrained <double,MyHS>::Functions::t ufns;
    ufns.f.reset(new Rosen);
    #if 0
    ustate.H_type=peopt::Operators::SR1;
    //ustate.H_type=peopt::Operators::External;
    peopt::Unconstrained <double,MyHS>::Functions
        ::init(peopt::Messaging(),ustate,ufns);
    #endif
    
    peopt::EqualityConstrained <double,MyHS,MyHS>::Functions::t efns;
    efns.f.reset(new Rosen);
    efns.g.reset(new Utility);
    estate.H_type=peopt::Operators::SR1;
    peopt::EqualityConstrained <double,MyHS,MyHS>::Functions
        ::init(peopt::Messaging(),estate,efns);
    
    peopt::InequalityConstrained <double,MyHS,MyHS>::Functions::t ifns;
    ifns.f.reset(new Rosen);
    ifns.h.reset(new Utility);
    #if 0
    istate.H_type=peopt::Operators::SR1;
    peopt::InequalityConstrained <double,MyHS,MyHS>::Functions
        ::init(peopt::Messaging(),istate,ifns);
    #endif
    
    peopt::Constrained <double,MyHS,MyHS,MyHS>::Functions::t cfns;
    cfns.f.reset(new Rosen);
    cfns.g.reset(new Utility);
    cfns.h.reset(new Utility);
    cstate.H_type=peopt::Operators::SR1;
    peopt::InequalityConstrained <double,MyHS,MyHS>::Functions
        ::init(peopt::Messaging(),cstate,cfns);
    
    // Construct the Rosenbrock fucntion
    Rosen f;

    // Do some finite difference tests on the Rosenbrock function
    peopt::Diagnostics::gradientCheck <> (peopt::Messaging(),f,x,dx);
    peopt::Diagnostics::hessianCheck <> (peopt::Messaging(),f,x,dx);
    peopt::Diagnostics::hessianSymmetryCheck <> (peopt::Messaging(),f,x,dx,dxx);
    
    // Construct the utility function
    Utility g;

    // Do some finite difference tests on the utility function
    peopt::Diagnostics::derivativeCheck <> (peopt::Messaging(),g,x,dx,dy);
    peopt::Diagnostics::derivativeAdjointCheck <>(peopt::Messaging(),g,x,dx,dy);
    peopt::Diagnostics::secondDerivativeCheck <> (peopt::Messaging(),g,x,dx,dy);

    // Setup the optimization problem
    #if 0
    // Newton's method
    ustate.H_type = peopt::Operators::External;
    ustate.iter_max = 100;
    ustate.eps_krylov = 1e-10;
    #endif

    // BFGS
    #if 0
    ustate.algorithm_class = peopt::AlgorithmClass::LineSearch;
    ustate.dir = peopt::LineSearchDirection::BFGS;
    ustate.stored_history = 10;
    ustate.iter_max = 300;
    #endif
    
    // Newton-CG 
    #if 1
    ustate.algorithm_class = peopt::AlgorithmClass::LineSearch;
    ustate.dir = peopt::LineSearchDirection::NewtonCG;
    ustate.H_type = peopt::Operators::External;
    ustate.eps_krylov = 1e-2;
    ustate.iter_max = 100;
    #endif

    // Solve the optimization problem
    peopt::Unconstrained <double,MyHS>::Algorithms
        ::getMin(peopt::Messaging(),ufns,ustate);

    // Print out the final answer
    const std::vector <double>& opt_x=*(ustate.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;
}
