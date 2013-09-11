// This example minimizes a simple quadratic function.  In order
// to modify the parameters, look at the file quadratic.param.  Only the
// Optizelle section is run, so copy in different parameter selections
// into that header in order to experiment.
// For reference, the optimal solution to this function is is (1,1,1).

#include <vector>
#include <iostream>
#include <string>
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define the quadratic function 
// 
// f(x,y,z)=(x-1)^2+(2y-2)^2+(3z-3)^2
//
struct Quad : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef Optizelle::Rm <double> X;

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
struct QuadHInv : public Optizelle::Operator <double,Optizelle::Rm,Optizelle::Rm> {
public:
    typedef Optizelle::Rm <double> X;
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

int main(int argc,char* argv[]){
    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "quadratic <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Generate an initial guess 
    std::vector <double> x(3);
    x[0]=-1.2; x[1]=1.1; x[2]=2.;

    // Create a direction for the finite difference tests
    std::vector <double> dx(3);
    dx[0]=-.5; dx[1]=.5; dx[2]=.5;
    
    // Create another direction for the finite difference tests
    std::vector <double> dxx(3);
    dxx[0]=.75; dxx[1]=.25;; dxx[2]=1.25;

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <double,Optizelle::Rm>::State::t state(x);

    // Read the parameters from file
    Optizelle::json::Unconstrained <double,Optizelle::Rm>::read(Optizelle::Messaging(),
        fname,state);

    // Create the bundle of functions 
    Optizelle::Unconstrained <double,Optizelle::Rm>::Functions::t fns;
    fns.f.reset(new Quad);
    fns.PH.reset(new QuadHInv(state.x.back())); 
    
    // Do some finite difference tests on the quadratic function
    Optizelle::Diagnostics::gradientCheck <> (Optizelle::Messaging(),*(fns.f),x,dx);
    Optizelle::Diagnostics::hessianCheck <> (Optizelle::Messaging(),*(fns.f),x,dx);
    Optizelle::Diagnostics::hessianSymmetryCheck <> (Optizelle::Messaging(),*(fns.f),
        x,dx,dxx);
    
    // Solve the optimization problem
    Optizelle::Unconstrained <double,Optizelle::Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Setup the optimization problem

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ',' << opt_x[2] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::Unconstrained <double,Optizelle::Rm>::write_restart(
        Optizelle::Messaging(),"quadratic.rst",state);

    // Successful termination
    return EXIT_SUCCESS;
}
