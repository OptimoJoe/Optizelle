// Optimize a simple problem with an optimal solution of (0.5,.25)

#include <iostream>
#include <iomanip>
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"

// Define a simple objective where 
// 
// f(x,y)=-x+y
//
struct MyObj : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef double Real;
    typedef Optizelle::Rm <Real> X;

    // Evaluation 
    double eval(const X::Vector& x) const {
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
// h(x,y) = [ y x ] >= 0
//          [ x 1 ]
//
struct MyIneq 
    : public Optizelle::VectorValuedFunction <double,Optizelle::Rm,Optizelle::SQL> 
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::SQL <double> Y;
    typedef double Real;

    // y=h(x) 
    void eval(
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y(1,1,1)=x[1];
        y(1,1,2)=x[0];
        y(1,2,1)=x[0];
        y(1,2,2)=Real(1.);
    }

    // y=h'(x)dx
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

    // z=h'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        //z[0]= Real(2.)*dy(1,1,2);
        z[0]= dy(1,1,2)+dy(1,2,1);
        z[1]= dy(1,1,1);
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
    // Create some type shortcuts
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::SQL <double> Z;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;

    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "simple_sdp_cone <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Generate an initial guess for the primal
    X_Vector x(2);
    x[0]=1.2; x[1]=3.1;

    // Generate an initial guess for the dual
    std::vector <Optizelle::Natural> sizes(1);
        sizes[0]=2;
    std::vector <Optizelle::Cone::t> types(1);
        types[0]=Optizelle::Cone::Semidefinite;
    Z_Vector z(types,sizes);
    Z::id(z);

    // Create an optimization state
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::State::t state(x,z);
    
    // Read the parameters from file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::read(Optizelle::Messaging(),fname,state);

    // Create a bundle of functions
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Solve the optimization problem
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::Algorithms::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << std::setprecision(16) << std::scientific 
        << "The optimal point is: (" << state.x[0] << ','
	<< state.x[1] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::write_restart(Optizelle::Messaging(),"solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
