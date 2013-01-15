#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// Test the Rosenbrock function with an optimal solution of (1,1)
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
    struct Rosen : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef peopt::Rm <double> X;

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

BOOST_AUTO_TEST_SUITE(rosenbrock)

BOOST_AUTO_TEST_CASE(newton_cg) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess for Rosenbrock 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::NewtonCG;
    state.H_type = peopt::Operators::UserDefined;
    state.eps_krylov = 1e-2;
    state.iter_max = 100;
    state.eps_grad = 1e-8;
    state.eps_dx = 1e-8;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,1), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 18 
    BOOST_CHECK(state.iter == 18);
}

BOOST_AUTO_TEST_CASE(tr_newton) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess for Rosenbrock 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 100;
    state.eps_krylov = 1e-10;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,1), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 21
    BOOST_CHECK(state.iter == 21);
}

BOOST_AUTO_TEST_CASE(bfgs) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess for Rosenbrock 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::BFGS;
    state.stored_history = 10;
    state.iter_max = 300;
    state.msg_level = 0;
    state.eps_dx = 1e-16;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,1), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 21
    BOOST_CHECK(state.iter == 24);
}

BOOST_AUTO_TEST_CASE(sr1) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess for Rosenbrock 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.H_type = peopt::Operators::SR1;
    state.stored_history = 5;  // Over sampling the Hessian makes a difference
    state.history_reset = 10;
    state.iter_max = 300;
    state.eps_krylov = 1e-16;
    state.eps_grad = 1e-10;
    state.eps_dx = 1e-10;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,1), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 51
    BOOST_CHECK(state.iter == 50);
}

BOOST_AUTO_TEST_SUITE_END()
