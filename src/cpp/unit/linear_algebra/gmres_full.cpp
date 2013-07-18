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
    Natural iter_max = 200;

    // Set how often we restart GMRES
    Natural rst_freq = 0;

    // Create some operator 
    BasicOperator <double> A(m);
    for(Natural i=1;i<=m*m;i++)
        A.A[i-1]=cos(pow(i,m-1));
    
    // Create some right hand side
    std::vector <double> b(m);
    for(Natural i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create the left preconditioner
    BasicOperator <double> Ml_inv(m);
    for(Natural i=1;i<=m*m;i++)
        Ml_inv.A[i-1]=cos(pow(30.+i,m-1));
    
    // Create the right preconditioner
    BasicOperator <double> Mr_inv(m);
    for(Natural i=1;i<=m*m;i++)
        Mr_inv.A[i-1]=cos(pow(55.+i,m-1));

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
    CHECK(err_iter.second == m);
    
    // Check the relative error between the true solution and that
    // returned from GMRES
    std::vector <double> x_star(5);
    x_star[0] = -1.203932331447497;
    x_star[1] = -0.186416740769010;
    x_star[2] = -0.457476984550115;
    x_star[3] = -0.830522778995837;
    x_star[4] = -1.125112777803922;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    CHECK(err < 1e-14);

    // Declare success
    return EXIT_SUCCESS;
}

