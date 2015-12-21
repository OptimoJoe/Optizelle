#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/linalg.h"
#include "linear_algebra.h"
#include "unit.h"

// In this problem, we have
// A = [ 1 -1 ]
//     [-1  1 ]
// b = [ 3 ]
//     [ 4 ]
// This has no solution.  On the first iteration, CG will move
// in the steepest descent direction, which is b.  In order to check the code
// for moving the center of a trust-region, we put the offset at [3;4] with 
// a radius of 7.5.  By setting the center in the opposite direction with a
// radius of 7.5, it should only move half the distance.
int main() {
    // Create a type shortcut
    typedef Optizelle::Rm <double> X;

    // Set the size of the problem
    Natural m = 2;

    // Set the stopping tolerance
    double eps_trunc = 1e-12;

    // Set the maximum number of iterations
    Natural iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 7.5;

    // Create some operator 
    BasicOperator <double> A(m);
    A.A[0]=1.;
    A.A[1]=-1.;
    A.A[2]=-1.;
    A.A[3]=1.;
    
    // Create some right hand side
    std::vector <double> b(m);
    b[0]=3.;
    b[1]=4.;
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;

    // Create a vector for the solution 
    std::vector <double> x(m);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Create a vector for the offset of the trust-region
    std::vector <double> x_offset(m);
    x_offset[0]=3.;
    x_offset[1]=4.;

    // Solve this linear system
    double residual_err0, residual_err; 
    Natural iter;
    Optizelle::TruncatedStop::t trunc_stop;
    auto safeguard_failed = Natural(0);
    auto alpha_safeguard = double(0.);
    Optizelle::truncated_cg <double,Optizelle::Rm>
        (A,b,W,eps_trunc,iter_max,1,delta,x_offset,false,1,
            no_safeguard <double,Optizelle::Rm>,x,x_cp,
            residual_err0,residual_err,iter,trunc_stop,safeguard_failed,
            alpha_safeguard);

    // Check that the size of x is 2.5 
    double norm_x = sqrt(X::innr(x,x));
    CHECK(std::abs(norm_x-2.5) < 1e-8);

    // Check that the solution is [1.5;2]
    std::vector <double> x_star(m);
    x_star[0] = 1.5; 
    x_star[1] = 2.; 
    std::vector <double> residual = x_star;
    Optizelle::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(Optizelle::Rm <double>::innr(residual,residual))
        /(1+sqrt(Optizelle::Rm <double>::innr(x_star,x_star)));
    CHECK(err < 1e-14);

    // Declare success
    return EXIT_SUCCESS;
}
