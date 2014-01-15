#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/linalg.h"
#include "linear_algebra.h"
#include "unit.h"

int main() {
    // Create a type shortcut
    typedef Optizelle::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    Natural m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    Natural iter_max = 200;

    // Set how often we restart GMRES
    Natural rst_freq = 0;

    // Create some operator with only three elements on the diagonal 
    BasicOperator <double> A(m);
    for(Natural i=1;i<=m*m;i++) A.A[i-1]=0;
    A.A[0]=2.;
    A.A[2+2*m]=3.;
    A.A[4+4*m]=4.;
    
    // Create some right hand side with ones in the places that correspond
    // to the operator above.  In this case, we have them at 1, 3, and 5.
    std::vector <double> b(m);
    for(Natural i=0;i<m;i++) b[i] = 0.;
    b[0]= 1.;
    b[2]= 1.;
    b[4]= 1.;
    
    // Create the left preconditioner.  Do this by inverting the matrix by
    // hand.
    BasicOperator <double> Ml_inv(m);
    for(Natural i=1;i<=m*m;i++) Ml_inv.A[i-1]=0.;
    Ml_inv.A[0]=1./2.;
    Ml_inv.A[2+2*m]=1./3.;
    Ml_inv.A[4+4*m]=1./4.;
    
    // Set the right preconditioner to the identity
    IdentityOperator <double> Mr_inv;

    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create an empty GMRES manipulator
    Optizelle::EmptyGMRESManipulator <double,Optizelle::Rm> gmanip;

    // Solve this linear system
    std::pair <double,Natural> err_iter = Optizelle::gmres <double,Optizelle::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    CHECK(err_iter.second == 1);
    
    // Check the relative error between the true solution and that
    // returned from GMRES
    std::vector <double> x_star(5);
    x_star[0] = 0.5;  
    x_star[1] = 0.; 
    x_star[2] = 1./3.; 
    x_star[3] = 0.; 
    x_star[4] = 0.25; 
    std::vector <double> residual = x_star;
    Optizelle::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(Optizelle::Rm <double>::innr(residual,residual))
        /(1+sqrt(Optizelle::Rm <double>::innr(x_star,x_star)));
    CHECK(err < 1e-14);
    
    // Declare success
    return EXIT_SUCCESS;
}

