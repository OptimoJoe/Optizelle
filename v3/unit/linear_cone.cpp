#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt.h"
#include "vspaces.h"

// Test a simple optimization problem with an optimal solution of (1/3,1/3)

namespace {
    // Squares its input
    template <typename Real>
    Real sq(Real x){
        return x*x; 
    }

    // Define a simple objective where 
    // 
    // f(x,y)=(x+1)^2+(y+1)^2
    //
    struct MyObj : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef peopt::Rm <double> X;
        typedef double Real;

        // Evaluation 
        double operator () (const X::Vector& x) const {
            return sq(x[0]+Real(1.))+sq(x[1]+Real(1.));
        }

        // Gradient
        void grad(
            const X::Vector& x,
            X::Vector& g
        ) const {
            g[0]=2*x[0]+2;
            g[1]=2*x[1]+2;
        }

        // Hessian-vector product
        void hessvec(
            const X::Vector& x,
            const X::Vector& dx,
            X::Vector& H_dx
        ) const {
            H_dx[0]= Real(2.)*dx[0]; 
            H_dx[1]= Real(2.)*dx[1]; 
        }
    };

    // Define simple inequalities 
    //
    // g(x,y)= [ x + 2y >= 1 ] 
    //         [ 2x + y >= 1 ] 
    //
    struct MyIneq
        : public peopt::VectorValuedFunction <double,peopt::Rm,peopt::Rm>
    {
        typedef peopt::Rm <double> X;
        typedef peopt::Rm <double> Y;
        typedef double Real;

        // y=f(x) 
        void operator () (
            const X::Vector& x,
            Y::Vector& y
        ) const {
            y[0]=x[0]+Real(2.)*x[1]-Real(1.);
            y[1]=Real(2.)*x[0]+x[1]-Real(1.);
        }

        // y=f'(x)dx
        void p(
            const X::Vector& x,
            const X::Vector& dx,
            Y::Vector& y
        ) const {
            y[0]= dx[0]+Real(2.)*dx[1];
            y[1]= Real(2.)*dx[0]+dx[1];
        }

        // z=f'(x)*dy
        void ps(
            const X::Vector& x,
            const Y::Vector& dy,
            X::Vector& z
        ) const {
            z[0]= dy[0]+Real(2.)*dy[1];
            z[1]= Real(2.)*dy[0]+dy[1];
        }

        // z=(f''(x)dx)*dy
        void pps(
            const X::Vector& x,
            const X::Vector& dx,
            const Y::Vector& dy,
            X::Vector& z
        ) const {
            X::zero(z);
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

BOOST_AUTO_TEST_SUITE(linear_cone)

BOOST_AUTO_TEST_CASE(newton_cg) {
    // Create a type shortcut
    using peopt::Rm;

    // Generate an initial guess for the primal
    std::vector <double> x(2);
    x[0]=2.1; x[1]=1.1;

    // Generate an initial guess for the dual
    std::vector <double> z(2);
    z[0]=1.; z[1]=1.;

    // Create an optimization state
    peopt::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);
    
    // Setup some parameters 
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::NewtonCG;
    state.H_type = peopt::Operators::External;
    state.eps_krylov = 1e-10;
    state.iter_max = 300;
    state.eps_s = 1e-16;
    state.eps_g = 1e-6;
    state.mu_trg=1e-6;
    state.sigma=0.05;
    state.gamma=.995;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);

    // Check the relative error between the true solution, (1/3,1/3), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1./3.; x_star[1]=1./3.;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 28);
}

BOOST_AUTO_TEST_CASE(tr_newton) {
    // Create a type shortcut
    using peopt::Rm;

    // Generate an initial guess for the primal
    std::vector <double> x(2);
    x[0]=2.1; x[1]=1.1;

    // Generate an initial guess for the dual
    std::vector <double> z(2);
    z[0]=1.; z[1]=1.;

    // Create an optimization state
    peopt::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);
    
    // Setup some parameters 
    state.H_type = peopt::Operators::External;
    state.iter_max = 100;
    state.eps_krylov = 1e-10;
    state.eps_s = 1e-16;
    state.eps_g = 1e-10;
    state.mu_trg=1e-6;
    state.sigma=0.05;
    state.gamma=.995;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);

    // Check the relative error between the true solution, (1/3,1/3), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1./3.; x_star[1]=1./3.;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 17);
}

BOOST_AUTO_TEST_CASE(bfgs) {
    // Create a type shortcut
    using peopt::Rm;

    // Generate an initial guess for the primal
    std::vector <double> x(2);
    x[0]=2.1; x[1]=1.1;

    // Generate an initial guess for the dual
    std::vector <double> z(2);
    z[0]=1.; z[1]=1.;

    // Create an optimization state
    peopt::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);
    
    // Setup some parameters 
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::BFGS;
    state.stored_history = 10;
    state.iter_max = 300;
    state.eps_s = 1e-16;
    state.mu_trg=1e-6;
    state.sigma=0.05;
    state.gamma=.995;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);

    // Check the relative error between the true solution, (1/3,1/3), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1./3.; x_star[1]=1./3.;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 54);
}

BOOST_AUTO_TEST_SUITE_END()
