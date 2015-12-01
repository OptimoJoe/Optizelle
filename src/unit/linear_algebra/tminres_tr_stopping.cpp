#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/linalg.h"
#include "linear_algebra.h"
#include "unit.h"

int main() {
    // Create a type shortcut
    typedef Optizelle::Rm <double> X;

    // Set the size of the problem
    Natural m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    Natural iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 0.1;

    // Create some operator 
    BasicOperator <double> A(m);
    for(Natural j=1;j<=m;j++)
        for(Natural i=1;i<=m;i++) {
            Natural I = j+(i-1)*m;
            Natural J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(Natural i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Create a vector for the center of the trust-region
    std::vector <double> x_cntr(m);
    Optizelle::Rm <double>::zero(x_cntr);

    // Solve this linear system
    double residual_err0, residual_err; 
    Natural iter;
    Optizelle::KrylovStop::t krylov_stop;
    auto failed_safeguard = Natural(0);
    auto alpha_safeguard = double(0.);
    Optizelle::truncated_minres <double,Optizelle::Rm>
        (A,b,W,eps_krylov,iter_max,1,delta,x_cntr,1,
            no_safeguard <double,Optizelle::Rm>,x,x_cp,
            residual_err0,residual_err,iter,krylov_stop,failed_safeguard,
            alpha_safeguard);

    // Check that the size of x is just the trust-region radius
    double norm_x = sqrt(X::innr(x,x));
    CHECK(std::abs(norm_x-delta) < 1e-8);

    // Declare success
    return EXIT_SUCCESS;
}
