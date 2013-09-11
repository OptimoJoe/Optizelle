// This tests our ability to do derivative and symmetry checks on scalar-valued
// functions.
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "unit.h"

namespace {
    // Squares its input
    template <typename Real>
    Real sq(Real x){
        return x*x; 
    }

    // Define the Rosenbrock function where
    // 
    // f(x,y)=(1-x)^2+100(y-x^2)^2
    //
    struct Rosen : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
        typedef Optizelle::Rm <double> X;

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
}

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

    // Do the finite difference tests
    double err=Optizelle::Diagnostics::gradientCheck <> (Optizelle::Messaging(),f,x,dx);
    CHECK(err < 1e-11);
    err=Optizelle::Diagnostics::hessianCheck <> (Optizelle::Messaging(),f,x,dx);
    CHECK(err < 1e-14);
    err=Optizelle::Diagnostics::hessianSymmetryCheck <>
        (Optizelle::Messaging(),f,x,dx,dxx);
    CHECK(err < 1e-12);

    // If we've made it this far, we're successful
    return EXIT_SUCCESS;

}
