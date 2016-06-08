// Solve the problem 
//
// min x + y + z
// st  [ x y; y z] >= 0 <==> xz >= y^2
//     y >=Q z <==> y >= |z|
//     z >= 1
//
// which should have an optimal solution of (1,1,1).  Basically, we have
// a problem with a semidefinite, quadratic, and linear cone.  Then, we're
// going to run a few iterations, stop the optimization, write the restart,
// load the restart, and then continue with the optimization.  This is basically
// a way to both show a bunch of features and add another unit test.

#include <iostream>
#include <iomanip>
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"

// Define a simple objective where 
//  
// f(x,y,z)=x+y+z
//
struct MyObj : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef double Real;
    typedef Optizelle::Rm <Real> X;
    
    // Evaluation 
    double eval(X::Vector const & x) const {
        auto sum = Real(0.);
        for(auto const & ele: x) 
            sum += ele;
        return sum;
    }

    // Gradient
    void grad(
        X::Vector const & x,
        X::Vector & grad 
    ) const {
        for(auto & ele : grad)
            ele = Real(1.);
    }

    // Hessian-vector product
    void hessvec(
        X::Vector const & x,
        X::Vector const & dx,
        X::Vector & H_dx
    ) const {
        X::zero(H_dx);
    }
};

// Define a simple SQL inequality 
//
// h(x,y,z) = [ x y ] >=S 0
//            [ y z ]
//            y >=Q z
//            z >=L 1
//
struct MyIneq 
    : public Optizelle::VectorValuedFunction
        <double,Optizelle::Rm,Optizelle::SQL> 
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::SQL <double> Z;
    typedef double Real;

    // z=h(x) 
    void eval(
        X::Vector const & x,
        Z::Vector & z
    ) const {
        z(1,1,1)=x[0];
        z(1,1,2)=x[1];
        z(1,2,1)=x[1];
        z(1,2,2)=x[2];

        z(2,1)=x[1];
        z(2,2)=x[2];

        z(3,1)=x[2]-1;
    }

    // y=h'(x)dx
    void p(
        X::Vector const & x,
        X::Vector const & dx,
        Z::Vector & z
    ) const {
        z(1,1,1)=dx[0];
        z(1,1,2)=dx[1];
        z(1,2,1)=dx[1];
        z(1,2,2)=dx[2];

        z(2,1)=dx[1];
        z(2,2)=dx[2];

        z(3,1)=dx[2];
    }

    // x_hat=h'(x)*dz
    void ps(
        X::Vector const & x,
        Y::Vector const & dz,
        X::Vector & x_hat 
    ) const {
        X::zero(x_hat);

        // Remember, the input to this function may not be symmetric, so
        // compute accordingly 
        x_hat[0] += dz(1,1,1);
        x_hat[1] += dz(1,1,2)+dz(1,2,1);
        x_hat[2] += dz(1,2,2);

        x_hat[1] += dz(2,1);
        x_hat[2] += dz(2,2);

        x_hat[2] += dz(3,1);
    }

    // x_hat=(h''(x)dx)*dy
    void pps(
        X::Vector const & x,
        X::Vector const & dx,
        Y::Vector const & dz,
        X::Vector & x_hat 
    ) const {
        X::zero(x_hat);
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
        std::cerr << "sql_restart <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    auto fname = argv[1];

    // Generate an initial guess for the primal
    auto x = std::vector <double> { 5.5, 3.3, 2.2};
    
    // Generate an initial guess for the dual
    auto z = Z_Vector(
        {Optizelle::Cone::Semidefinite,
         Optizelle::Cone::Quadratic,
         Optizelle::Cone::Linear},
        {2,2,1});

    // Create an optimization state
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::State::t state(x,z);
    
    // Read the parameters from file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::read(fname,state);

    // Create a bundle of functions
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.h.reset(new MyIneq);

    // Solve the optimization problem
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::Algorithms::getMin(Optizelle::Messaging::stdout,fns,state);
    
    // Write an intermediate restart file 
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::write_restart("restart.json",state);

    // Read in the restart file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::read_restart("restart.json",x,z,state);

    // Change the maximum number of iterations to something larger and reset
    // the convergence flag
    state.iter_max = 500;
    state.opt_stop = Optizelle::OptimizationStop::NotConverged;

    // Finish solving the optimization problem
    Optizelle::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::Algorithms::getMin(Optizelle::Messaging::stdout,fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::OptimizationStop::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << std::setprecision(16) << std::scientific 
        << "The optimal point is: (" << state.x[0] << ','
        << state.x[1] << ',' << state.x[2] << ')' << std::endl;

    // Write out the final answer to file
    Optizelle::json::InequalityConstrained <double,Optizelle::Rm,Optizelle::SQL>
        ::write_restart("solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
