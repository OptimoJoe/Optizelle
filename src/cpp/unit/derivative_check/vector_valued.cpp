#include "peopt/peopt.h"
#include "peopt/vspaces.h"
#include "unit.h"

// This tests our ability to do derivative and symmetry checks on vector-valued
// functions

namespace {
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

    // Define some utility function where
    //
    // g(x,y)= [ cos(x)sin(y) ]
    //         [ 3 x^2 y + y^3]
    //         [ log(x) + 3y^5]
    //
    struct Utility  : public peopt::VectorValuedFunction
        <double,peopt::Rm,peopt::Rm>
    {
        typedef peopt::Rm <double> X;
        typedef peopt::Rm <double> Y;

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
}

int main() {
    // Create some arbitrary vectors 
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    std::vector <double> dx(2);
    dx[0] = 3.4; dx[1] = 4.3;
    std::vector <double> dy(3);
    dy[0]=.3; dy[1]=-.12; dy[2]=1.2;

    // Construct the utility function
    Utility g;

    // Do the finite difference tests
    double err=peopt::Diagnostics::derivativeCheck <>
        (peopt::Messaging(),g,x,dx,dy);
    CHECK(err < 1e-12);
    err=peopt::Diagnostics::derivativeAdjointCheck <>
        (peopt::Messaging(),g,x,dx,dy);
    CHECK(err < 1e-15);
    err=peopt::Diagnostics::secondDerivativeCheck <>
        (peopt::Messaging(),g,x,dx,dy);
    CHECK(err < 1e-12);

    // If we've made it this far, we're successful
    return EXIT_SUCCESS;
}
