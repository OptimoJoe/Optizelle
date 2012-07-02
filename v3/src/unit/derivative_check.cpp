#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// This tests our ability to do derivative and symmetry checks

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

    // Define the Rosenbrock function where
    // 
    // f(x,y)=(1-x)^2+100(y-x^2)^2
    //
    struct Rosen : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef peopt::Rm <double> X;

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
    struct Utility  : public peopt::VectorValuedFunction
        <double,peopt::Rm,peopt::Rm>
    {
        typedef peopt::Rm <double> X;
        typedef peopt::Rm <double> Y;

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

    // Does not output anything to the user unless its an error 
    struct SilentMessaging : public peopt::Messaging {
        // Prints a message
        void print(const std::string msg) const { }

        // Prints an error
        void error(const std::string msg) const {
            std::cerr << msg << std::endl;
            exit(EXIT_FAILURE);
        }
    };
}

BOOST_AUTO_TEST_SUITE(derivative_check)

BOOST_AUTO_TEST_CASE(scalar_valued_function) {
    // Create some arbitrary vectors in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    std::vector <double> dx(2);
    dx[0] = 3.4; dx[1] = 4.3;
    std::vector <double> dxx(2);
    dxx[0] = 3.2; dxx[1] = 1.1;

    // Construct the Rosenbrock function
    Rosen f;

    // Do the finite difference tests
    double err=peopt::Diagnostics::gradientCheck <> (SilentMessaging(),f,x,dx);
    BOOST_CHECK(err < 1e-11);
    err=peopt::Diagnostics::hessianCheck <> (SilentMessaging(),f,x,dx);
    BOOST_CHECK(err < 1e-14);
    err=peopt::Diagnostics::hessianSymmetryCheck <>
        (SilentMessaging(),f,x,dx,dxx);
    BOOST_CHECK(err < 1e-12);

}

BOOST_AUTO_TEST_CASE(vector_valued_function) {
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
        (SilentMessaging(),g,x,dx,dy);
    BOOST_CHECK(err < 1e-12);
    err=peopt::Diagnostics::derivativeAdjointCheck <>
        (SilentMessaging(),g,x,dx,dy);
    BOOST_CHECK(err < 1e-15);
    err=peopt::Diagnostics::secondDerivativeCheck <>
        (SilentMessaging(),g,x,dx,dy);
    BOOST_CHECK(err < 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
