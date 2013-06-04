// This example minimizes the Rosenbrock function using the sample
// vector spaces and parsing routines provided by peopt.  In order
// to modify the parameters, look at the file run.peopt.  Only the
// peopt section is run, so copy in different parameter selections
// into that header in order to experiment.
// For reference, the optimal solution to the Rosenbrock function
// is (1,1).

#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"
#include "peopt/json.h"

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosen
    : public peopt::ScalarValuedFunction <double,peopt::Rm>
{
    typedef peopt::Rm <double> X;

    // Evaluation of the Rosenbrock function
    double operator () (const X::Vector& x) const {
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]));
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& grad
    ) const {
        grad[0]=-400.*x[0]*(x[1]-sq(x[0]))-2.*(1.-x[0]);
        grad[1]=200.*(x[1]-sq(x[0]));
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
    	H_dx[0]=(1200.*sq(x[0])-400.*x[1]+2)*dx[0]-400.*x[0]*dx[1];
        H_dx[1]=-400.*x[0]*dx[0]+200.*dx[1];
    }
};

// Define a perfect preconditioner for the Hessian
struct RosenHInv : public peopt::Operator <double,peopt::Rm,peopt::Rm> {
public:
    typedef peopt::Rm <double> X;
    typedef X::Vector X_Vector;
private:
    X_Vector& x;
public:
    RosenHInv(X::Vector& x_) : x(x_) {}
    void operator () (const X_Vector& dx,X_Vector &result) const {
        double one_over_det=1./(400000.*x[0]*x[0]-80000.*x[1]+400.);
        result[0]=one_over_det*(200.*dx[0]+400.*x[0]*dx[1]);
        result[1]=one_over_det*
            (400.*x[0]*dx[0]+(1200.*x[0]*x[0]-400.*x[1]+2.)*dx[1]);
    }
};

int main(){
    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=-1.2; x[1]=1.;

    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,peopt::Rm>::State::t state(x);

    // Read the parameters from file
    peopt::json::Unconstrained <double,peopt::Rm>
        ::read(peopt::Messaging(),"preconditioned_nonlinear_cg.peopt",state);

    // Create the bundle of functions 
    peopt::Unconstrained <double,peopt::Rm>::Functions::t fns;
    fns.f.reset(new Rosen);
    fns.PH.reset(new RosenHInv(state.x.back()));
    
    // Solve the optimization problem
    peopt::Unconstrained <double,peopt::Rm>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        peopt::StoppingCondition::to_string(state.opt_stop) <<
        std::endl;

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;

    // Successful termination
    return EXIT_SUCCESS;
}
