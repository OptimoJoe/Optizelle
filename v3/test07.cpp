#include <vector>
#include <iostream>
#include <string>
#include "peopt.h"
#include "vspaces.h"
#include "json.h"

#define MyHS peopt::Rm

// Squares its input
template <typename Real>
Real sq(Real x){
    return x*x; 
}

// Define the Rosenbrock function where
// 
// f(x,y)=(1-x)^2+100(y-x^2)^2
//
struct Rosen : public peopt::ScalarValuedFunction <double,MyHS> {
    typedef MyHS <double> X;

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

int main(int argc,char* argv[]){

    // Parse our inputs
    bool output_params=false;
    std::string fname;
    bool usage_problem=false;;
    if(argc >=1 && argc <=3) {
        // Loop over all the arguments
        for(int i=1;i<argc;i++) {
            // Get the argument
            std::string arg=argv[i];

            // Check if we're outputing our parameters
            if(arg=="-params") output_params=true;

            // Otherwise, we have a filename.  
            else {
                // Save the filename.  If we already have a filename, there's
                // a problem.
                if(fname.size()==0)
                    fname=arg;
                else
                    usage_problem=true;

                // Check if the file is good
                std::ifstream test(arg.c_str());
                if(!test.good()) usage_problem=true; 
                test.close();
            }
        }
    }
    if(usage_problem) {
        std::cout << "Usage: test07 [-params] [param_file]" << std::endl;
        return EXIT_FAILURE;
    } 

    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=-1.2; x[1]=1.;

    // Create a direction for the finite difference tests
    std::vector <double> dx(2);
    dx[0]=-.5; dx[1]=.5;
    
    // Create another direction for the finite difference tests
    std::vector <double> dxx(2);
    dxx[0]=.75; dxx[1]=.25;

    // Create an optimization state
    peopt::Unconstrained <double,MyHS>::State::t state(x);

    // If we have parameters, read in from file
    if(fname.size() > 0) 
        peopt::json::Unconstrained <double,MyHS>::read(peopt::Messaging(),
            fname,state);

    // Print parameters to screen if requested
    if(output_params) {
        std::cout << peopt::json::Unconstrained <double,MyHS>::to_string(state);
        return EXIT_SUCCESS;
    }

    // Create a bundle of functions
    peopt::Unconstrained <double,MyHS>::Functions::t fns;
    fns.f.reset(new Rosen);
    
    // Do some finite difference tests on the Rosenbrock function
    peopt::Diagnostics::gradientCheck <> (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianCheck <> (peopt::Messaging(),*(fns.f),x,dx);
    peopt::Diagnostics::hessianSymmetryCheck <> (peopt::Messaging(),*(fns.f),x,
        dx,dxx);
    
    // Solve the optimization problem
    peopt::Unconstrained <double,MyHS>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        peopt::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.x.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
	<< opt_x[1] << ')' << std::endl;

    return EXIT_SUCCESS;
}
