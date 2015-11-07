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
    double delta = 100.;

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

    // Get the norm of the RHS
    double norm_b = std::sqrt(X::innr(b,b));
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;

    // Create an initial guess at the solution
    std::vector <double> x(m);
    for(Natural i=1;i<=m;i++) x[i-1]=1.;

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Create a vector for the center of the trust-region
    std::vector <double> x_cntr(m);
    Optizelle::Rm <double>::zero(x_cntr);

    // Solve this linear system
    double residual_err0, residual_err; 
    Natural iter;
    Optizelle::KrylovStop::t krylov_stop;
    Optizelle::truncated_cd <double,Optizelle::Rm>
        (A,b,W,eps_krylov,iter_max,1,delta,x_cntr,true,x,x_cp,
            residual_err0,residual_err,iter,krylov_stop);

    // Check the error is less than our tolerance 
    CHECK(residual_err < eps_krylov*norm_b);

    // Check that we ran to the maximum number of iterations
    CHECK(iter == m);
    
    // Check the relative error between the true solution and that
    // returned from TPCG 
    std::vector <double> x_star(5);
    x_star[0] = 0.062210523692158425;
    x_star[1] = -0.027548098303754341;
    x_star[2] = -0.11729291808469694;
    x_star[3] = -0.080812473373141375;
    x_star[4] = 0.032637688404329734;
    std::vector <double> residual = x_star;
    Optizelle::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(Optizelle::Rm <double>::innr(residual,residual))
        /(1+sqrt(Optizelle::Rm <double>::innr(x_star,x_star)));
    CHECK(err < 1e-14);

    // Check that the returned solution is different than the Cauchy point
    Optizelle::Rm <double>::copy(x_cp,residual);
    Optizelle::Rm <double>::axpy(-1,x,residual);
    err=std::sqrt(Optizelle::Rm <double>::innr(residual,residual))
        /(1+sqrt(Optizelle::Rm <double>::innr(x_cp,x_cp)));
    CHECK(err > 1e-4);

    // Declare success
    return EXIT_SUCCESS;
}
