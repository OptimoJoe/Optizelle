// A problem that helps us determine how to scale inequality constrained
// optimization problems

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

// Grab Optizelle's natural type
typedef Optizelle::Natural Natural;

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Define a simple objective where 
// 
// f(x) = 0.5 || x - c ||^2
//
struct MyObj : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef Optizelle::Rm <double> X;

    // Center of the quadratic objective
    X::Vector const & c;

    // Grab the center during constuction
    MyObj(X::Vector const & c_) : c(c_) {}

    // Evaluation 
    double eval(const X::Vector& x) const {
        double acc(0.);
        for(Natural i=0;i<x.size();i++)
            acc+=sq(x[i]-c[i]);
        return 0.5*acc;
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& grad
    ) const {
        for(Natural i=0;i<x.size();i++)
            grad[i]=x[i]-c[i];
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
        for(Natural i=0;i<x.size();i++)
            H_dx[i]=dx[i];
    }
};

// Define simple inequalities 
//
// h(x) = x - lb 
//
struct MyIneq
    :public Optizelle::VectorValuedFunction<double,Optizelle::Rm,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;

    // Lower bound 
    X::Vector const & lb;
    
    // Grab the lower bound during constuction
    MyIneq(X::Vector const & lb_) : lb(lb_) {}

    // y=h(x) 
    void eval(
        const X::Vector& x,
        Y::Vector& y
    ) const {
        for(Natural i=0;i<x.size();i++)
            y[i]=x[i]-lb[i];
    }

    // y=h'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        for(Natural i=0;i<x.size();i++)
            y[i]=dx[i];
    }

    // z=h'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        for(Natural i=0;i<x.size();i++)
            z[i]=dy[i];
    }

    // z=(h''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        X::zero(z);
    }
};

int main(int argc,char* argv[]){
    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "inequality_scaling <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Create a type shortcut
    using Optizelle::Rm;

    // Generate an initial guess
    std::vector <double> x(10);
    for(Natural i=0;i<x.size();i++)
        x[i]=1. + pow(.1,i);

    // Allocate memory for the inequality multipler 
    std::vector <double> z(10);

    // Create the center of the objective function
    std::vector <double> c(10);
    for(Natural i=0;i<c.size();i++)
        c[i]=-1.;

    // Create the lower bound for the problem
    std::vector <double> lb(10);
    for(Natural i=0;i<lb.size();i++)
        lb[i]=1.;

    // Create an optimization state
    Optizelle::InequalityConstrained <double,Rm,Rm>::State::t
        state(x,z);

    // Read the parameters from file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::read(Optizelle::Messaging(),fname,state);
    
    // Create a bundle of functions
    Optizelle::InequalityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj(c));
    fns.h.reset(new MyIneq(lb));

    // Solve the optimization problem
    Optizelle::InequalityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    std::cout << "The optimal point is: [" << std::endl;
    for(Natural i=0;i<state.x.size();i++)
        std::cout << std::scientific << std::setprecision(16)
            << state.x[i] << std::endl;
    std::cout << "]" << std::endl;

    // Write out the final answer to file
    Optizelle::json::InequalityConstrained<double,Rm,Rm>
        ::write_restart(Optizelle::Messaging(),"solution.json",state);

    // Return that the program exited properly
    return EXIT_SUCCESS;
}
