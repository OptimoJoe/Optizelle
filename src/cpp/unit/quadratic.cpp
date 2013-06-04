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
    // f(x,y,z)=(x-1)^2+(2y-2)^2+(3z-3)^2
    //
    // This has a minimum at (1,1,1).
    struct Quad : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef peopt::Rm <double> X;

        // Evaluation of the quadratic function
        double operator () (const X::Vector& x) const {
            return sq(x[0]-1.)+sq(2*x[1]-2.)+sq(3*x[2]-3.);
        }

        // Gradient
        void grad(
            const X::Vector& x,
            X::Vector& g
        ) const {
            g[0]=2*x[0]-2;
            g[1]=8*x[1]-8;
            g[2]=18*x[2]-18;
        }

        // Hessian-vector product
        void hessvec(
            const X::Vector& x,
            const X::Vector& dx,
            X::Vector& H_dx
        ) const {
            H_dx[0]= 2*dx[0];
            H_dx[1]= 8*dx[1]; 
            H_dx[2]= 18*dx[2]; 
        }
    };

    // Define an almost perfect preconditioner for the Hessian
    struct QuadHInv : public peopt::Operator <double,peopt::Rm,peopt::Rm> {
    public:
        typedef peopt::Rm <double> X;
        typedef X::Vector X_Vector;
    private:
        X_Vector& x;
    public:
        QuadHInv(X::Vector& x_) : x(x_) {}
        void operator () (const X_Vector& dx,X_Vector &result) const {
            result[0]=dx[0]/2.;
            result[1]=dx[1]/8.;
            result[2]=dx[2];
        }
    };
}

BOOST_AUTO_TEST_SUITE(quadratic)

BOOST_AUTO_TEST_CASE(tr_newton) {
    // Create a type shortcut
    using peopt::Rm;

    // Generate an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.eps_krylov = 1e-10;
    state.delta = 100.;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  Since this is quadratic, Newton's
    // method should work in a single iteration.
    BOOST_CHECK(state.iter == 2);
}

BOOST_AUTO_TEST_CASE(ncg_fletcher_reeves) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::FletcherReeves;
    state.linesearch_iter_max=200;
    state.iter_max = 50;
    state.eps_grad = 1e-7;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  With an exact line-search, this should
    // be the number of variables+1.
    BOOST_CHECK(state.iter == 4);
}

BOOST_AUTO_TEST_CASE(ncg_polak_ribiere){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::PolakRibiere;
    state.iter_max = 50;
    state.linesearch_iter_max=200;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  With an exact line-search, this should
    // be the number of variables+1.
    BOOST_CHECK(state.iter == 4);
}

BOOST_AUTO_TEST_CASE(ncg_hestenes_stiefel){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::HestenesStiefel;
    state.iter_max = 50;
    state.linesearch_iter_max=200;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  With an exact line-search, this should
    // be the number of variables+1.
    BOOST_CHECK(state.iter == 4);
}

BOOST_AUTO_TEST_CASE(precon_ncg_fletcher_reeves) {
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::FletcherReeves;
    state.PH_type = peopt::Operators::UserDefined;
    state.linesearch_iter_max=200;
    state.iter_max = 50;
    state.eps_grad = 1e-7;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);
    fns.PH.reset(new QuadHInv(state.x.back())); 

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  We knock out two eigenvalues with
    // the preconditioner, so this should be the number of (variables+1)-1.
    BOOST_CHECK(state.iter == 3);
}

BOOST_AUTO_TEST_CASE(precon_ncg_polak_ribiere){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::PolakRibiere;
    state.PH_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.linesearch_iter_max=200;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);
    fns.PH.reset(new QuadHInv(state.x.back())); 

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  We knock out two eigenvalues with
    // the preconditioner, so this should be the number of (variables+1)-1.
    BOOST_CHECK(state.iter == 3);
}

BOOST_AUTO_TEST_CASE(precon_ncg_hestenes_stiefel){
    // Create a type shortcut
    using peopt::Rm;

    // Create an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;
    
    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Setup some algorithmic parameters
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::HestenesStiefel;
    state.PH_type = peopt::Operators::UserDefined;
    state.iter_max = 50;
    state.linesearch_iter_max=200;
    state.msg_level = 0;
    
    // Create the bundle of functions 
    peopt::Unconstrained <double,Rm>::Functions::t fns;
    fns.f.reset(new Quad);
    fns.PH.reset(new QuadHInv(state.x.back())); 

    // Solve the optimization problem
    peopt::Unconstrained <double,Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);
    
    // Check the relative error between the true solution and that found in 
    // the state
    std::vector <double> x_star(3);
    x_star[0] = 1.0; x_star[1]=1.0; x_star[2]=1.0;
    std::vector <double> residual = x_star;
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations.  We knock out two eigenvalues with
    // the preconditioner, so this should be the number of (variables+1)-1.
    BOOST_CHECK(state.iter == 3);
}

BOOST_AUTO_TEST_SUITE_END()
