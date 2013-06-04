#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// Define a simple optimization problem with an optimal solution of (0.5,0.25)
namespace {
    // Squares its input
    template <typename Real>
    Real sq(Real x){
        return x*x; 
    }

    // Define a simple objective where 
    // 
    // f(x,y)=-x+y
    //
    struct MyObj : public peopt::ScalarValuedFunction <double,peopt::Rm> {
        typedef double Real;
        typedef peopt::Rm <Real> X;

        // Evaluation 
        double operator () (const X::Vector& x) const {
            return -x[0]+x[1]; 
        }

        // Gradient
        void grad(
            const X::Vector& x,
            X::Vector& g
        ) const {
            g[0]=Real(-1.);
            g[1]=Real(1.);
        }

        // Hessian-vector product
        void hessvec(
            const X::Vector& x,
            const X::Vector& dx,
            X::Vector& H_dx
        ) const {
            H_dx[0]= Real(0.);
            H_dx[1]= Real(0.);
        }
    };

    // Define a simple SDP inequality 
    //
    // g(x,y) = [ y x ] >= 0
    //          [ x 1 ]
    //
    struct MyIneq 
        : public peopt::VectorValuedFunction <double,peopt::Rm,peopt::SQL> 
    {
        typedef peopt::Rm <double> X;
        typedef peopt::SQL <double> Y;
        typedef double Real;

        // y=f(x) 
        void operator () (
            const X::Vector& x,
            Y::Vector& y
        ) const {
            y(1,1,1)=x[1];
            y(1,1,2)=x[0];
            y(1,2,1)=x[0];
            y(1,2,2)=Real(1.);
        }

        // y=f'(x)dx
        void p(
            const X::Vector& x,
            const X::Vector& dx,
            Y::Vector& y
        ) const {
            y(1,1,1)=dx[1];
            y(1,1,2)=dx[0];
            y(1,2,1)=dx[0];
            y(1,2,2)=Real(0.);
        }

        // z=f'(x)*dy
        void ps(
            const X::Vector& x,
            const Y::Vector& dy,
            X::Vector& z
        ) const {
            z[0]= dy(1,1,2)+dy(1,2,1);
            z[1]= dy(1,1,1);
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
}


BOOST_AUTO_TEST_SUITE(sdp_cone)

BOOST_AUTO_TEST_CASE(newton_cg) {
    
    // Create some type shortcuts
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <peopt::Natural> sizes(1); sizes[0]=2;
    std::vector <peopt::Cone::t> types(1); types[0]=peopt::Cone::Semidefinite;
    Z_Vector z(peopt::Messaging(),types,sizes);
    Z::id(z);

    // Create an optimization state
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Setup the optimization problem
    state.algorithm_class = peopt::AlgorithmClass::LineSearch;
    state.dir = peopt::LineSearchDirection::NewtonCG;
    state.H_type = peopt::Operators::UserDefined;
    state.eps_krylov = 1e-10;
    state.iter_max = 100;
    state.msg_level = 0;
    state.eps_dx = 1e-16;
    state.eps_grad = 1e-13;
    state.eps_mu = 1e-11;
    state.sigma = 0.1;
    state.gamma = 0.99;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Check the relative error between the true solution, (0.5,0.25), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 0.5; x_star[1]=0.25;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 15);
}

BOOST_AUTO_TEST_CASE(tr_newton) {
    
    // Create some type shortcuts
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <peopt::Natural> sizes(1); sizes[0]=2;
    std::vector <peopt::Cone::t> types(1); types[0]=peopt::Cone::Semidefinite;
    Z_Vector z(peopt::Messaging(),types,sizes);
    Z::id(z);

    // Create an optimization state
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Setup the optimization problem
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 100;
    state.msg_level = 0;
    state.eps_krylov = 1e-10;
    state.eps_dx = 1e-16;
    state.eps_grad = 1e-10;
    state.eps_mu = 1e-11;
    state.sigma = 0.1;
    state.gamma = 0.99;
    state.delta = 100;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Check the relative error between the true solution, (0.5,0.25), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 0.5; x_star[1]=0.25;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 13);
}

BOOST_AUTO_TEST_CASE(tr_newton_predictor_corrector) {
    
    // Create some type shortcuts
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <peopt::Natural> sizes(1); sizes[0]=2;
    std::vector <peopt::Cone::t> types(1); types[0]=peopt::Cone::Semidefinite;
    Z_Vector z(peopt::Messaging(),types,sizes);
    Z::id(z);

    // Create an optimization state
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Setup the optimization problem
    state.H_type = peopt::Operators::UserDefined;
    state.iter_max = 100;
    state.msg_level = 0;
    state.eps_krylov = 1e-10;
    state.eps_dx = 1e-16;
    state.eps_grad = 1e-8;
    state.eps_mu = 1e-8;
    state.gamma = 0.9;
    state.sigma = 0.5;
    state.delta = 100;
    state.cstrat = peopt::CentralityStrategy::PredictorCorrector;

    // Solve the optimization problem
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Check the relative error between the true solution, (0.5,0.25), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 0.5; x_star[1]=0.25;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 17);
}


BOOST_AUTO_TEST_CASE(sr1) {
    
    // Create some type shortcuts
    typedef peopt::Rm <double> X;
    typedef peopt::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <peopt::Natural> sizes(1); sizes[0]=2;
    std::vector <peopt::Cone::t> types(1); types[0]=peopt::Cone::Semidefinite;
    Z_Vector z(peopt::Messaging(),types,sizes);

    // Create an optimization state
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Create a bundle of functions
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Setup the optimization problem
    state.algorithm_class = peopt::AlgorithmClass::TrustRegion;
    state.H_type = peopt::Operators::SR1;
    state.stored_history = 2;
    state.eps_krylov=1e-10;
    state.iter_max = 100;
    state.msg_level = 0;
    state.sigma = 0.5;
    state.gamma = 0.90;
    state.eps_dx = 1e-16;
    state.eps_grad = 1e-9;
    state.eps_mu = 1e-10;
    state.delta = 100;
    
    // Solve the optimization problem
    peopt::InequalityConstrained <double,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Check the relative error between the true solution, (0.5,0.25), and that
    // found in the state
    std::vector <double> x_star(2);
    x_star[0] = 0.5; x_star[1]=0.25;
    std::vector <double> residual = x_star;
    peopt::Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(peopt::Rm <double>::innr(residual,residual))
        /(1+sqrt(peopt::Rm <double>::innr(x_star,x_star)));
    BOOST_CHECK(err < 1e-6);

    // Check the number of iterations 
    BOOST_CHECK(state.iter == 29);
}


BOOST_AUTO_TEST_SUITE_END()
