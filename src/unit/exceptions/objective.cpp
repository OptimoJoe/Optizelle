//  Tests our ability to throw a native exception and catch it 

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"

// Define a function that just throws an exception
struct Objective:public Optizelle::ScalarValuedFunction <double,Optizelle::Rm>{
    typedef Optizelle::Rm <double> X;

    // Evaluation
    double eval(const X::Vector & x) const {
        throw std::runtime_error("Objective"); 
    }

    // Gradient
    void grad(
        X::Vector const & x,
        X::Vector & grad
    ) const {
        throw std::runtime_error("Gradient"); 
    }

    // Hessian-vector product
    void hessvec(
        X::Vector const & x,
        X::Vector const & dx,
        X::Vector & H_dx
    ) const {
        throw std::runtime_error("Hessian-vector product"); 
    }
};

int main() try {
    // Create a type shortcut
    using Optizelle::Rm;

    // Allocate memory for an initial guess
    auto x = std::vector <double> {1.2, 2.3};
    
    // Create an optimization state
    Optizelle::Unconstrained<double,Rm>::State::t state(x);

    // Create a bundle of functions
    Optizelle::Unconstrained<double,Rm>::Functions::t fns;
    fns.f.reset(new Objective);

    // Try to catch the error
    Optizelle::Unconstrained <double,Rm>::Algorithms
        ::getMin(Optizelle::Messaging::stdout,fns,state);

} catch(std::runtime_error const & e) {
    if(std::string(e.what()) != "Gradient")
        throw;
}
