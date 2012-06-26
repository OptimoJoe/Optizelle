#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt.h"
#include "vspaces.h"

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
    state.H_type = peopt::Operators::External;
    state.eps_krylov = 1e-2;
    state.iter_max = 100;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);
    
    // Check the relative error between the true solution, (1,1), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 22
    BOOST_CHECK(state.iter == 22);
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
    state.H_type = peopt::Operators::External;
    state.iter_max = 100;
    state.eps_krylov = 1e-10;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);
    
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
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Rosen);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(SilentMessaging(),fns,state);
    
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

BOOST_AUTO_TEST_SUITE_END()
