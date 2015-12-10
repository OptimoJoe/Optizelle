// This example minimizes a simple quadratic function.  For reference, the 
// optimal solution to this function is (1,1,1).

#include <vector>
#include <iostream>
#include <string>
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Define the quadratic function 
// 
// f(x,y,z)=(x-1)^2+(2y-2)^2+(3z-3)^2
//
struct Quad : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef Optizelle::Rm <double> X;

    // Evaluation of the quadratic function
    double eval(const X::Vector& x) const {
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
    void eval(const X_Vector& dx,X_Vector &result) const {
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

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <double,Optizelle::Rm>::State::t state(x);

    // Read the parameters from file
    Optizelle::json::Unconstrained <double,Optizelle::Rm>::read(
        Optizelle::Messaging(),fname,state);

    // Create the bundle of functions 
    Optizelle::Unconstrained <double,Optizelle::Rm>::Functions::t fns;
    fns.f.reset(new Quad);
    fns.PH.reset(new QuadHInv(state.x)); 

    // Solve the optimization problem
    Optizelle::Unconstrained <double,Optizelle::Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ',' << state.x[2] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::Unconstrained <double,Optizelle::Rm>::write_restart(
        Optizelle::Messaging(),"solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
