#include <vector>
#include <iostream>
#include <string>
#include "peopt.h"

// Defines the vector space used for optimization.
struct MyOps : public peopt::Operations <std::vector <double>, double> {
    typedef std::vector <double> Vector;
    typedef double Real;

    // Create an empty, uninitialized vector
    Vector create() const {
        return std::vector <double> ();
    }
    
    // Memory allocation and size setting
    void init(const Vector& x, Vector& y) const {
        y.resize(x.size());
    }

    // y <- x (Shallow.  No memory allocation.)
    void copy(const Vector& x, Vector& y) const {
        for(unsigned int i=0;i<x.size();i++){
            y[i]=x[i];
        }
    }

    // x <- alpha * x
    void scal(const Real& alpha, Vector& x) const {
        for(unsigned int i=0;i<x.size();i++){
            x[i]=alpha*x[i];
        }
    }

    // x <- 0 
    void zero(Vector& x) const {
        for(unsigned int i=0;i<x.size();i++){
            x[i]=0.;
        }
    }

    // y <- alpha * x + y
    void axpy(const Real& alpha, const Vector& x, Vector& y) const {
        for(unsigned int i=0;i<x.size();i++){
            y[i]=alpha*x[i]+y[i];
        }
    }

    // innr <- <x,y>
    Real innr(const Vector& x,const Vector& y) const {
        Real z=0;
        for(unsigned int i=0;i<x.size();i++){
            z+=x[i]*y[i];
        }
        return z;
    }
};

// Squares its input
double sq(double x){
    return x*x; 
}

// Cubes its input
double cub(double x){
    return x*x*x; 
}

// Quads its input
double quad(double x){
    return x*x*x*x; 
}

// Quints its input
double quint(double x){
    return x*x*x*x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosen : public peopt::ScalarValuedFunction <std::vector<double>,double> {
    // Evaluation of the Rosenbrock function
    double operator () (const std::vector<double>& x) const {
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]));
    }

    // Gradient
    void grad(const std::vector<double>& x,std::vector<double>& g) const {
        g[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0]);
        g[1]=200*(x[1]-sq(x[0]));
    }

    // Hessian-vector product
    void hessvec(
        const std::vector<double>& x,
        const std::vector<double>& dx,
        std::vector<double>& H_dx 
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
struct Utility  : public peopt::VectorValuedFunction
    <std::vector<double>,std::vector<double>,double> {
    typedef std::vector<double> Domain;
    typedef std::vector<double> Codomain;
    typedef double Real;

    // y=f(x) 
    void operator () (const Domain& x,Codomain& y) const {
        y[0]=cos(x[0])*sin(x[1]);
        y[1]=Real(3)*sq(x[0])*x[1]+cub(x[1]);
        y[2]=log(x[0])+Real(3.)*quint(x[1]);
    }

    // y=f'(x)dx
    void p(const Domain& x,const Domain& dx,Codomain& y) const {
        y[0]= -sin(x[0])*sin(x[1])*dx[0]
              +cos(x[0])*cos(x[1])*dx[1];
        y[1]= Real(6.)*x[0]*x[1]*dx[0]
              +(Real(3.)*sq(x[0])+Real(3.)*sq(x[1]))*dx[1];
        y[2]= Real(1.)/x[0]*dx[0]
              +Real(15.)*quad(x[1])*dx[1];
    }

    // z=f'(x)*dy
    void ps(const Domain& x,const Codomain& dy,Domain& z) const {
        z[0]= -sin(x[0])*sin(x[1])*dy[0]
              +Real(6.)*x[0]*x[1]*dy[1]
              +Real(1.)/x[0]*dy[2];
        z[1]= cos(x[0])*cos(x[1])*dy[0]
              +(Real(3.)*sq(x[0])+Real(3.)*sq(x[1]))*dy[1]
              +Real(15.)*quad(x[1])*dy[2];
    }
    // z=(f''(x)dx)*dy
    void pps(const Domain& x,const Domain& dx,
        const Codomain& dy,Domain& z) const {
        z[0] = (-cos(x[0])*dx[0]*sin(x[1])-sin(x[0])*cos(x[1])*dx[1])*dy[0]
               +(Real(6.)*dx[0]*x[1] + Real(6.)*x[0]*dx[1])*dy[1]
               +(Real(-1.)/sq(x[0])*dx[0])*dy[2];
        z[1] = (-sin(x[0])*dx[0]*cos(x[1])-cos(x[0])*sin(x[1])*dx[1])*dy[0]
               +(Real(6.)*x[0]*dx[0]+Real(6.)*x[1]*dx[1])*dy[1]
               +(Real(60.)*cub(x[1])*dx[1])*dy[2];
    }
};


int main(){

    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=1.3; x[1]=1.;

    // Create a direction for the finite difference tests
    std::vector <double> dx(2);
    dx[0]=-.5; dx[1]=.5;
    
    // Create another direction for the finite difference tests
    std::vector <double> dxx(2);
    dxx[0]=.75; dxx[1]=.25;

    // Create a vector in the codomain of the utility function
    std::vector <double> dy(3);
    dy[0]=.3; dy[1]=-.12; dy[2]=1.2;

    // Construct the Rosenbrock fucntion
    Rosen f;

    // Do some finite difference tests on the Rosenbrock function
    peopt::gradientCheck <> (peopt::Messaging(),MyOps(),f,x,dx);
    peopt::hessianCheck <> (peopt::Messaging(),MyOps(),f,x,dx);
    peopt::hessianSymmetryCheck <> (peopt::Messaging(),MyOps(),f,x,dx,dxx);
    
    // Construct the utility function
    Utility g;

    // Do some finite difference tests on the utility function
    peopt::derivativeCheck <> (peopt::Messaging(),MyOps(),MyOps(),g,x,dx,dy);
    peopt::derivativeAdjointCheck <> (peopt::Messaging(),MyOps(),MyOps(),
        g,x,dx,dy);
    peopt::secondDerivativeCheck <> (peopt::Messaging(),MyOps(),MyOps(),
        g,x,dx,dy);

  
#if 0
    // Create a state and setup the problem
    peopt::core<MyVS>::State state(x);


    // General problem setup
    state.eps_g=1e-10;
    state.eps_s=1e-10;
    state.iter_max=200;
    state.eps_krylov=1e-8;

    // Newton's method
    state.H_type=peopt::Operators::External;
    state.algorithm_class=peopt::AlgorithmClass::TrustRegion;

    // BFGS
    //state.dir=peopt::LineSearchDirection::BFGS;
    //state.kind=peopt::LineSearchKind::GoldenSection;
    //state.algorithm_class=peopt::AlgorithmClass::LineSearch;
    //state.stored_history=10;

    // Create a function, gradient, and Hessian for this problem
    RosenObjective F;
    RosenGradient G;
    // Make sure to link the Hessian to the current optimization
    // iterate
    RosenHessian H(*(state.u.begin()));
   
    // Do a finite difference test for the gradient and Hessian
    std::cout << "Gradient finite difference check" << std::endl;
    peopt::derivativeCheck<MyVS>(F,G,x,eta);

    std::cout << "Hessian finite difference check" << std::endl;
    peopt::derivativeCheck<MyVS,MyVS>(G,H,x,eta,x);

    // Optimize the problem
    peopt::core<MyVS>::getMin(state,F,G,H);

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.u.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;
#endif
}
