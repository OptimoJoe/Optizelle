#include "peopt/peopt.h"
#include "peopt/vspaces.h"
#include "peopt/linalg.h"
#include "linear_algebra.h"
#include "unit.h"

int main() {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    Natural m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    Natural iter_max = 300;

    // Set how often we restart GMRES
    Natural rst_freq = 3;

    // Create some operator 
    BasicOperator <double> A(m);
    for(Natural i=1;i<=m*m;i++)
        A.A[i-1]=cos(pow(i,2));
    
    // Create some right hand side
    std::vector <double> b(m);
    for(Natural i=1;i<=m;i++) b[i-1] = cos(i+25); 

    // Create the left preconditioner
    IdentityOperator <double> Ml_inv;
    
    // Create the right preconditioner
    IdentityOperator <double> Mr_inv;

    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create an empty GMRES manipulator
    peopt::GMRESManipulator <double,peopt::Rm> gmanip;

    // Solve this linear system
    std::pair <double,Natural> err_iter = peopt::gmres <double,peopt::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    CHECK(err_iter.second == 242);
    
    // Declare success
    return EXIT_SUCCESS;
}

