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
    
    // Create someright hand side
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
    
    // Create someright hand side with ones in the places that correspond
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
    
    // Create someright hand side with ones in the places that correspond
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
    
    // Create someright hand side
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



BOOST_AUTO_TEST_SUITE_END()
