// This example demonstrates how to run a series of diagnostic tests
// on functions and then immediately exit.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"

// Grab Optizelle's Natural type
using Optizelle::Natural;

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
struct Rosenbrock :
    public Optizelle::ScalarValuedFunction <double,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;

    // Evaluation of the Rosenbrock function
    double eval(const X::Vector & x) const {
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

// Define some utility function where
//
// g(x)= [ cos(x1) sin(x2)   ]
//       [ 3 x1^2 x2 + x2 ^3 ]
//       [ log(x1) + 3 x2 ^5 ]
//
struct Utility  : public Optizelle::VectorValuedFunction
    <double,Optizelle::Rm,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;

    // y=g(x) 
    void eval(
        const X::Vector & x,
        Y::Vector & y
    ) const {
        y[0]=cos(x[0])*sin(x[1]);
        y[1]=3.*sq(x[0])*x[1]+cub(x[1]);
        y[2]=log(x[0])+3.*quint(x[1]);
    }

    // y=g'(x)dx
    void p(
        const X::Vector & x,
        const X::Vector & dx,
        Y::Vector & y
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
        const X::Vector & x,
        const Y::Vector & dy,
        X::Vector & z
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
        const X::Vector & x,
        const X::Vector & dx,
        const Y::Vector & dy,
        X::Vector & z
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

    // Create a type shortcut
    using Optizelle::Rm;

    // Allocate memory for an initial guess and equality multiplier 
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    std::vector <double> z(3);
    
    // Create an optimization state
    Optizelle::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Modify the state so that we just run our diagnostics and exit
    state.dscheme = Optizelle::DiagnosticScheme::DiagnosticsOnly;
    state.f_diag = Optizelle::FunctionDiagnostics::SecondOrder;
    state.x_diag = Optizelle::VectorSpaceDiagnostics::Basic;
    state.h_diag = Optizelle::FunctionDiagnostics::SecondOrder;
    state.z_diag = Optizelle::VectorSpaceDiagnostics::EuclideanJordan;
    state.L_diag = Optizelle::FunctionDiagnostics::SecondOrder;
    
    // Create a bundle of functions
    Optizelle::InequalityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new Rosenbrock);
    fns.h.reset(new Utility);

    // Even though this looks like we're solving an optimization problem,
    // we're actually just going to run our diagnostics and then exit.
    Optizelle::InequalityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);
}
