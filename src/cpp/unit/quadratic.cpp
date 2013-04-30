#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// Test the behavior of different algorithms on an easy quadratic
namespace {
    // Squares its input
    template <typename Real>
    Real sq(Real x){
        return x*x; 
    }

    // Define the quadratic function 
    // 
    // f(x,y)=(x-1)^2+(x-2)^2
    //
    struct Quad : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef peopt::Rm <double> X;

        // Evaluation of the quadratic function 
        double operator () (const X::Vector& x) const {
            return sq(x[0]-1.)+sq(x[1]-2.);
        }

        // Gradient
        void grad(
            const X::Vector& x,
            X::Vector& g
        ) const {
            g[0]=2*x[0]-2;
            g[1]=2*x[1]-4;
        }

        // Hessian-vector product
        void hessvec(
            const X::Vector& x,
            const X::Vector& dx,
            X::Vector& H_dx
        ) const {
            H_dx[0]= 2*dx[0];
            H_dx[1]= 2*dx[1]; 
        }
    };
}

BOOST_AUTO_TEST_SUITE(quadratic)

BOOST_AUTO_TEST_CASE(tr_newton) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.eps_krylov = 1e-10;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,2), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=2.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check that the number of iterations is 1
    BOOST_CHECK(state.iter == 3);
}

BOOST_AUTO_TEST_CASE(ncg_fletcher_reeves) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::FletcherReeves;
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,2), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=2.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 5);
}

BOOST_AUTO_TEST_CASE(ncg_polak_ribiere){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::PolakRibiere;
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,2), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=2.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 6);
}

BOOST_AUTO_TEST_CASE(ncg_hestenes_stiefel){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(2);
    x[0] = -1.2; x[1] = 1.; 
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::HestenesStiefel;
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution, (1,2), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 1.0; x_star[1]=2.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 9);

    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state0(x);

    // Setup some algorithmic parameters
    state0.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state0.dir = peopt::LineSearchDirection::HestenesStiefel;
    state0.H_type = peopt::Operators::UserDefined;
    state0.iter_max = 50;
    state0.msg_level = 0;
    state0.linesearch_iter_max = 10;

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state0);

    // Check the number of iterations.  This should be lower than above.
    BOOST_CHECK(state0.iter == 6);
}

BOOST_AUTO_TEST_SUITE_END()
