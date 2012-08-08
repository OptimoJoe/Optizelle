#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"
#include "peopt/linalg.h"

// This tests our linear algebra routines 

template <typename Real>
struct BasicOperator : public peopt::Operator <Real,peopt::Rm,peopt::Rm> {
private:
    // Create some type shortcuts
    typedef std::vector <Real> X_Vector;
    typedef std::vector <Real> Y_Vector;
    
    // Size of the matrix
    unsigned int m;
public:
    // Storage for the matrix
    std::vector <Real> A;

    // Create an empty matrix, which must be filled by the user 
    BasicOperator(const unsigned int m_) : m(m_) {
        A.resize(m*m); 
    }
    
    // Apply the random matrix to the vector 
    void operator () (const X_Vector& x,Y_Vector &y) const  {
        for(unsigned int i=0;i<m;i++){
            y[i]=0;
            for(unsigned int j=0;j<m;j++)
                y[i]+=A[i+m*j]*x[j];
        }
    }
};

// Create the identity operator
template <typename Real>
struct IdentityOperator : public peopt::Operator <Real,peopt::Rm,peopt::Rm> {
private:
    // Create some type shortcuts
    typedef std::vector <Real> X_Vector;
    typedef std::vector <Real> Y_Vector;

public:
    // We don't need to do anything on creation
    IdentityOperator() {};

    // Just copy the input to the output 
    void operator () (const X_Vector& x,Y_Vector &y) const  {
        peopt::Rm <Real>::copy(x,y);
    }
};


BOOST_AUTO_TEST_SUITE(linear_algebra)

BOOST_AUTO_TEST_CASE(gmres_full) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set how often we restart GMRES
    unsigned int rst_freq = 0;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int i=1;i<=m*m;i++)
        A.A[i-1]=cos(pow(i,m-1));
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create the left preconditioner
    BasicOperator <double> Ml_inv(m);
    for(int i=1;i<=m*m;i++)
        Ml_inv.A[i-1]=cos(pow(30.+i,m-1));
    
    // Create the right preconditioner
    BasicOperator <double> Mr_inv(m);
    for(int i=1;i<=m*m;i++)
        Mr_inv.A[i-1]=cos(pow(55.+i,m-1));

    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create an empty GMRES manipulator
    peopt::GMRESManipulator <double,peopt::Rm> gmanip;

    // Solve this linear system
    std::pair <double,unsigned int> err_iter = peopt::gmres <double,peopt::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == m);
    
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
    BOOST_CHECK(err < 1e-14);
}

BOOST_AUTO_TEST_CASE(gmres_left_preconditioner) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set how often we restart GMRES
    unsigned int rst_freq = 0;

    // Create some operator with only three elements on the diagonal 
    BasicOperator <double> A(m);
    for(int i=1;i<=m*m;i++) A.A[i-1]=0;
    A.A[0]=2.;
    A.A[2+2*m]=3.;
    A.A[4+4*m]=4.;
    
    // Create some right hand side with ones in the places that correspond
    // to the operator above.  In this case, we have them at 1, 3, and 5.
    std::vector <double> b(m);
    for(unsigned int i=0;i<m;i++) b[i] = 0.;
    b[0]= 1.;
    b[2]= 1.;
    b[4]= 1.;
    
    // Create the left preconditioner.  Do this by inverting the matrix by
    // hand.
    BasicOperator <double> Ml_inv(m);
    for(int i=1;i<=m*m;i++) Ml_inv.A[i-1]=0.;
    Ml_inv.A[0]=1./2.;
    Ml_inv.A[2+2*m]=1./3.;
    Ml_inv.A[4+4*m]=1./4.;
    
    // Set the right preconditioner to the identity
    IdentityOperator <double> Mr_inv;

    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create an empty GMRES manipulator
    peopt::GMRESManipulator <double,peopt::Rm> gmanip;

    // Solve this linear system
    std::pair <double,unsigned int> err_iter = peopt::gmres <double,peopt::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == 1);
    
    // Check the relative error between the true solution and that
    // returned from GMRES
    std::vector <double> x_star(5);
    x_star[0] = 0.5;  
    x_star[1] = 0.; 
    x_star[2] = 1./3.; 
    x_star[3] = 0.; 
    x_star[4] = 0.25; 
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);
}

BOOST_AUTO_TEST_CASE(gmres_right_preconditioner) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set how often we restart GMRES
    unsigned int rst_freq = 0;

    // Create some operator with only three elements on the diagonal 
    BasicOperator <double> A(m);
    for(int i=1;i<=m*m;i++) A.A[i-1]=0;
    A.A[0]=2.;
    A.A[2+2*m]=3.;
    A.A[4+4*m]=4.;
    
    // Create some right hand side with ones in the places that correspond
    // to the operator above.  In this case, we have them at 1, 3, and 5.
    std::vector <double> b(m);
    for(unsigned int i=0;i<m;i++) b[i] = 0.;
    b[0]= 1.;
    b[2]= 1.;
    b[4]= 1.;
    
    // Create the right preconditioner.  Do this by inverting the matrix by
    // hand.
    BasicOperator <double> Mr_inv(m);
    for(int i=1;i<=m*m;i++) Mr_inv.A[i-1]=0.;
    Mr_inv.A[0]=1./2.;
    Mr_inv.A[2+2*m]=1./3.;
    Mr_inv.A[4+4*m]=1./4.;
    
    // Set the left preconditioner to the identity
    IdentityOperator <double> Ml_inv;

    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create an empty GMRES manipulator
    peopt::GMRESManipulator <double,peopt::Rm> gmanip;

    // Solve this linear system
    std::pair <double,unsigned int> err_iter = peopt::gmres <double,peopt::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == 1);
    
    // Check the relative error between the true solution and that
    // returned from GMRES
    std::vector <double> x_star(5);
    x_star[0] = 0.5;  
    x_star[1] = 0.; 
    x_star[2] = 1./3.; 
    x_star[3] = 0.; 
    x_star[4] = 0.25; 
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);
}

BOOST_AUTO_TEST_CASE(gmres_restart) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 300;

    // Set how often we restart GMRES
    unsigned int rst_freq = 3;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int i=1;i<=m*m;i++)
        A.A[i-1]=cos(pow(i,2));
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
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
    std::pair <double,unsigned int> err_iter = peopt::gmres <double,peopt::Rm>
        (A,b,eps_krylov,iter_max,rst_freq,Ml_inv,Mr_inv,gmanip,x);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == 242);
}

BOOST_AUTO_TEST_CASE(tpcd_basic_solve) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 100.;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == m);
    
    // Check the relative error between the true solution and that
    // returned from TPCG 
    std::vector <double> x_star(5);
    x_star[0] = 0.062210523692158425;
    x_star[1] = -0.027548098303754341;
    x_star[2] = -0.11729291808469694;
    x_star[3] = -0.080812473373141375;
    x_star[4] = 0.032637688404329734;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);

    // Check that the returned solution is different than the Cauchy point
    peopt::Rm <double>::copy(x_cp,residual);
    peopt::Rm <double>::axpy(-1,x,residual);
    err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_cp,x_cp)));
    BOOST_CHECK(err > 1e-4);
}

BOOST_AUTO_TEST_CASE(tpcd_tr_stopping) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 0.1;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check that the size of x is just the trust-region radius
    double norm_x = sqrt(X::innr(x,x));
    BOOST_CHECK_CLOSE(norm_x,delta,1e-8);
}

BOOST_AUTO_TEST_CASE(tpcd_cp) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 1;

    // Set the trust-reregion radius 
    double delta = 100.;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == 1);
    
    // Check that the returned solution and the Cauchy point are the same
    std::vector <double> residual = x_cp;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_cp,x_cp)));
    BOOST_CHECK(err < 1e-14);
}

BOOST_AUTO_TEST_CASE(tpcd_nullspace_solve) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 100.;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create a simple nullspace projector.  This projects out the first
    // two elements
    BasicOperator <double> W(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            W.A[I-1]=(i==j && i<=2) ? 1. : 0.;
        }

    
    // Create some right hand side.  Make sure that this is in the range
    // of A*W.
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = A.A[i-1]+A.A[i-1+m];
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    X::zero (x);

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we completed in two iterations.  This is due to the
    // nullspace projection
    BOOST_CHECK(err_iter.second == 2);
    
    // Check the relative error between the true solution and that
    // returned from TPCG.  Technically, these may not match exactly,
    // but Wx should equal x_star.
    std::vector <double> x_star(5);
    x_star[0] = 1.0; 
    x_star[1] = 1.0; 
    x_star[2] = 0.0;
    x_star[3] = 0.0;
    x_star[4] = 0.0;
    std::vector <double> residual = x_star;
    std::vector <double> Wx(m);
    W(x,Wx);
    peopt::Rm <double>::axpy(-1,Wx,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);

    // Check that the returned solution is different than the Cauchy point
    peopt::Rm <double>::copy(x_cp,residual);
    peopt::Rm <double>::axpy(-1,x,residual);
    err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_cp,x_cp)));
    BOOST_CHECK(err > 1e-4);
}

BOOST_AUTO_TEST_CASE(tpcd_starting_solution) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 100.;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    for(int i=1;i<=m;i++) x[i-1]=1.;

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // Check that we ran to the maximum number of iterations
    BOOST_CHECK(err_iter.second == m);
    
    // Check the relative error between the true solution and that
    // returned from TPCG 
    std::vector <double> x_star(5);
    x_star[0] = 0.062210523692158425;
    x_star[1] = -0.027548098303754341;
    x_star[2] = -0.11729291808469694;
    x_star[3] = -0.080812473373141375;
    x_star[4] = 0.032637688404329734;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);

    // Check that the returned solution is different than the Cauchy point
    peopt::Rm <double>::copy(x_cp,residual);
    peopt::Rm <double>::axpy(-1,x,residual);
    err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_cp,x_cp)));
    BOOST_CHECK(err > 1e-4);
}

BOOST_AUTO_TEST_CASE(tpcd_optimal_starting_solution) {
    // Create a type shortcut
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;

    // Set the size of the problem
    unsigned int m = 5;

    // Set the stopping tolerance
    double eps_krylov = 1e-12;

    // Set the maximum number of iterations
    unsigned int iter_max = 200;

    // Set the trust-reregion radius 
    double delta = 100.;

    // Create some operator 
    BasicOperator <double> A(m);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=m;i++) {
            unsigned I = j+(i-1)*m;
            unsigned J = i+(j-1)*m;
            if(i>j) {
                A.A[I-1]=cos(pow(I,m-1));
                A.A[J-1]=A.A[I-1];
            } else if(i==j)
                A.A[I-1]=cos(pow(I,m-1))+10;
        }
    
    // Create some right hand side
    std::vector <double> b(m);
    for(unsigned int i=1;i<=m;i++) b[i-1] = cos(i+25); 
    
    // Create some empty null-space projection 
    IdentityOperator <double> W;
    
    // Create an initial guess at the solution
    std::vector <double> x(m);
    x[0] = 0.062210523692158425;
    x[1] = -0.027548098303754341;
    x[2] = -0.11729291808469694;
    x[3] = -0.080812473373141375;
    x[4] = 0.032637688404329734;

    // Create a vector for the Cauchy point
    std::vector <double> x_cp(m);

    // Solve this linear system
    std::pair <double,unsigned int> err_iter
        = peopt::truncated_pcd <double,peopt::Rm>
            (A,b,W,eps_krylov,iter_max,delta,x,x_cp);

    // Check the error is less than our tolerance 
    BOOST_CHECK(err_iter.first < eps_krylov);

    // This one is tricky.  Ideally, we shouldn't do any work, but since
    // we only check for relative improvement in the residual, we'll actually
    // do full iterations since the residual isn't going to drop.
    BOOST_CHECK(err_iter.second == m);
    
    // Check the relative error between the true solution and that
    // returned from TPCG 
    std::vector <double> x_star(5);
    x_star[0] = 0.062210523692158425;
    x_star[1] = -0.027548098303754341;
    x_star[2] = -0.11729291808469694;
    x_star[3] = -0.080812473373141375;
    x_star[4] = 0.032637688404329734;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,x,residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-14);

    // Check that the returned solution and the Cauchy point are the same
    residual = x_cp;
    peopt::Rm <double>::axpy(-1,x,residual);
    err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_cp,x_cp)));
    BOOST_CHECK(err < 1e-14);
}

BOOST_AUTO_TEST_SUITE_END()
