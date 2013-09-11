// This example demonstrates how to create a generic vector space,
// define a scalar valued function (for the objective) and a vector valued
// function (for the constraints), and then run a series of diagnostic
// checks on them.
#include "optizelle/optizelle.h"

// Defines the vector space used for optimization.
template <typename Real>
struct MyHS { 
    typedef std::vector <Real> Vector;
    
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
struct Rosen : public Optizelle::ScalarValuedFunction <double,MyHS> {
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
// g(x)= [ cos(x1) sin(x2)   ]
//       [ 3 x1^2 x2 + x2 ^3 ]
//       [ log(x1) + 3 x2 ^5 ]
//
struct Utility  : public Optizelle::VectorValuedFunction
    <double,MyHS,MyHS>
{
    typedef MyHS <double> X;
    typedef MyHS <double> Y;

    // y=g(x) 
    void operator () (
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[0]=cos(x[0])*sin(x[1]);
        y[1]=3.*sq(x[0])*x[1]+cub(x[1]);
        y[2]=log(x[0])+3.*quint(x[1]);
    }

    // y=g'(x)dx
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

    // z=g'(x)*dy
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

    // z=(g''(x)dx)*dy
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

int main() {

    // Create some arbitrary vectors in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    std::vector <double> dx(2);
    dx[0] = 3.4; dx[1] = 4.3;
    std::vector <double> dxx(2);
    dxx[0] = 3.2; dxx[1] = 1.1;

    // Construct the Rosenbrock function
    Rosen f;

    // Do two finite difference checks and then check the symmetry of the
    // Hessian
    Optizelle::Diagnostics::gradientCheck <> (Optizelle::Messaging(),f,x,dx);
    Optizelle::Diagnostics::hessianCheck <> (Optizelle::Messaging(),f,x,dx);
    Optizelle::Diagnostics::hessianSymmetryCheck <> (Optizelle::Messaging(),f,x,dx,dxx);

    // Create some vectors in R^3 for testing the vector-valued function. 
    std::vector <double> dy(3);
    dy[0]=.3; dy[1]=-.12; dy[2]=1.2;

    // Construct the utility function
    Utility g;

    // Do the finite difference tests and check whether or not the
    // derivative is adjoint to the derivative adjoint.
    Optizelle::Diagnostics::derivativeCheck <> (Optizelle::Messaging(),g,x,dx,dy);
    Optizelle::Diagnostics::derivativeAdjointCheck <>(Optizelle::Messaging(),g,x,dx,dy);
    Optizelle::Diagnostics::secondDerivativeCheck <>(Optizelle::Messaging(),g,x,dx,dy);
}
