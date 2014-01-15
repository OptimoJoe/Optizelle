// Exercise 18.2 in Numerical Optimization by Nocedal and Wright.  This has
// an optimal solution of x = (-1.71,1.59,1.82.-0.763,-0.763).

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include <iostream>
#include <iomanip>

// Grab the Natural type from Optizelle
using Optizelle::Natural;

// Squares its input
template <typename Real>
Real sq(Real const & x){
    return x*x; 
}

// Indexing for vectors
Natural itok(Natural const & i) {
    return i-1;
}

// Indexing for packed storage
Natural ijtokp(Natural i,Natural j) {
    if(i>j) {
        Natural tmp=i;
        i=j;
        j=tmp;
    }
    return (i-1)+j*(j-1)/2;
}
    
// Indexing function for dense matrices 
Natural ijtok(Natural const & i,Natural const & j,Natural const & m) {
    return (i-1)+(j-1)*m;
}

// Indexing function for dense tensors 
Natural ijktol(Natural const & i,Natural const & j,Natural const & k,
    Natural const & m,Natural const & n
) {
    return (i-1)+(j-1)*m+(k-1)*m*n;
}

// 
// f(x) = exp(x1 x2 x3 x4 x5) - (1/2) (x1^3 + x2^3 + 1)^2
//
struct MyObj : public Optizelle::ScalarValuedFunction <double,Optizelle::Rm> {
    typedef Optizelle::Rm <double> X;
    typedef double Real;

    // Evaluation 
    double eval(const X::Vector& x) const {
        return exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - sq(pow(x[itok(1)],3)+pow(x[itok(2)],3)+Real(1.))/Real(2.);
    }

    // Gradient
    void grad(
        const X::Vector& x,
        X::Vector& g
    ) const {
        g[itok(1)]= x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - Real(3.)*sq(x[itok(1)])
            *(pow(x[itok(1)],3) + pow(x[itok(2)],3)+Real(1.));
        g[itok(2)]= x[itok(1)]*x[itok(3)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - Real(3.)*sq(x[itok(2)])
            * (pow(x[itok(1)],3) + pow(x[itok(2)],3) + Real(1.));
        g[itok(3)]= x[itok(1)]*x[itok(2)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);
        g[itok(4)]= x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);
        g[itok(5)] = x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);
    }

    // Hessian-vector product
    void hessvec(
        const X::Vector& x,
        const X::Vector& dx,
        X::Vector& H_dx
    ) const {
        // Allocate memory for the dense Hessian in packed storage
        std::vector <Real> H(15);

        // Compute the dense Hessian
        H[ijtokp(1,1)] =
            sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                - Real(6.)*x[itok(1)]*(pow(x[itok(1)],3)
                + pow(x[itok(2)],3)+Real(1.))
                - Real(9.)*pow(x[itok(1)],4);
        H[ijtokp(1,2)]=
            x[itok(1)]*x[itok(2)]*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(3)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)])
                - Real(9.)*sq(x[itok(1)])*sq(x[itok(2)]);
        H[ijtokp(1,3)]=
            x[itok(1)]*sq(x[itok(2)])*x[itok(3)]*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(1,4)]=
            x[itok(1)]*sq(x[itok(2)])*sq(x[itok(3)])*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(3)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(1,5)]=
            x[itok(1)]*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(3)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(2,2)] =
            sq(x[itok(1)])*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                - Real(6.)*x[itok(2)]
                * (pow(x[itok(1)],3)+pow(x[itok(2)],3)+Real(1.))
                - Real(9.)*pow(x[itok(2)],4);
        H[ijtokp(2,3)]=
            sq(x[itok(1)])*x[itok(2)]*x[itok(3)]*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(2,4)]=
            sq(x[itok(1)])*x[itok(2)]*sq(x[itok(3)])*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(3)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(2,5)]=
            sq(x[itok(1)])*x[itok(2)]*sq(x[itok(3)])*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(3)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(3,3)]=
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(4)])*sq(x[itok(5)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(3,4)]=
            sq(x[itok(1)])*sq(x[itok(2)])*x[itok(3)]*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(3,5)]=
            sq(x[itok(1)])*sq(x[itok(2)])*x[itok(3)]*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(4,4)]=
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(5)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(4,5)]=
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*x[itok(4)]*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(3)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]);
        H[ijtokp(5,5)]=
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]);

        // Compute the Hessian-vector product
        X::zero(H_dx);
        for(Natural i=1;i<=5;i++) {
            for(Natural j=1;j<=5;j++) {
                H_dx[i-1] += H[ijtokp(i,j)]*dx[j-1];
            }
        }
    }
};

//
// g(x)= [ x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 ]
//       [ x2 x3 - 5 x4 x5                       ]
//       [ x1^3 + x2^3 + 1                       ]
//
struct MyEq
    : public Optizelle::VectorValuedFunction <double,Optizelle::Rm,Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;
    typedef double Real;

    // y=g(x) 
    void eval(
        const X::Vector& x,
        Y::Vector& y
    ) const {
        y[itok(1)] = sq(x[itok(1)]) + sq(x[itok(2)]) + sq(x[itok(3)])
            + sq(x[itok(4)]) + sq(x[itok(5)]) - Real(10.);
        y[itok(2)] = x[itok(2)]*x[itok(3)] - Real(5.)*x[itok(4)]*x[itok(5)];
        y[itok(3)] = pow(x[itok(1)],3) + pow(x[itok(2)],3) + Real(1.);
    }

    // Generate a dense version of the Jacobian
    static void generateJac(const std::vector<Real>& x,std::vector <Real>& jac){
        jac[ijtok(1,1,3)] = Real(2.)*x[itok(1)];
        jac[ijtok(1,2,3)] = Real(2.)*x[itok(2)];
        jac[ijtok(1,3,3)] = Real(2.)*x[itok(3)];
        jac[ijtok(1,4,3)] = Real(2.)*x[itok(4)];
        jac[ijtok(1,5,3)] = Real(2.)*x[itok(5)];
        
        jac[ijtok(2,2,3)] = x[itok(3)];
        jac[ijtok(2,3,3)] = x[itok(2)];
        jac[ijtok(2,4,3)] = Real(-5.)*x[itok(5)];
        jac[ijtok(2,5,3)] = Real(-5.)*x[itok(4)];
        
        jac[ijtok(3,1,3)] = Real(3.)*sq(x[itok(1)]);
        jac[ijtok(3,2,3)] = Real(3.)*sq(x[itok(2)]);
    }

    // y=g'(x)dx
    void p(
        const X::Vector& x,
        const X::Vector& dx,
        Y::Vector& y
    ) const {
        // Generate a dense matrix that holds the Jacobian
        std::vector <Real> jac(15,Real(0.));

        // Compute a dense form of the Jacobian
        generateJac(x,jac);

        // Compute the Jacobian-vector product
        X::zero(y);
        for(Natural i=1;i<=3;i++) {
            for(Natural j=1;j<=5;j++) {
                y[itok(i)] += jac[ijtok(i,j,3)]*dx[itok(j)];
            }
        }
    }

    // z=g'(x)*dy
    void ps(
        const X::Vector& x,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        // Generate a dense matrix that holds the Jacobian
        std::vector <Real> jac(15,Real(0.));

        // Compute a dense form of the Jacobian
        generateJac(x,jac);

        // Compute the Jacobian transpose-vector product
        X::zero(z);
        for(Natural i=1;i<=3;i++) {
            for(Natural j=1;j<=5;j++) {
                z[itok(j)] += jac[ijtok(i,j,3)]*dy[itok(i)];
            }
        }
    }

    // z=(g''(x)dx)*dy
    void pps(
        const X::Vector& x,
        const X::Vector& dx,
        const Y::Vector& dy,
        X::Vector& z
    ) const {
        // Generate a dense tensor that holds the second derivative adjoNatural 
        std::vector <Real> D(75,Real(0.));
        D[ijktol(1,1,1,3,5)] = Real(2.);
        D[ijktol(1,2,2,3,5)] = Real(2.);
        D[ijktol(1,3,3,3,5)] = Real(2.);
        D[ijktol(1,4,4,3,5)] = Real(2.);
        D[ijktol(1,5,5,3,5)] = Real(2.);
        
        D[ijktol(2,2,3,3,5)] = Real(1.);
        D[ijktol(2,3,2,3,5)] = Real(1.);
        D[ijktol(2,4,5,3,5)] = Real(-5.);
        D[ijktol(2,5,4,3,5)] = Real(-5.);

        D[ijktol(3,1,1,3,5)] = Real(6.)*x[itok(1)];
        D[ijktol(3,2,2,3,5)] = Real(6.)*x[itok(2)];

        // Compute the action of this operator on our directions
        X::zero(z);
        for(Natural i=1;i<=3;i++) {
            for(Natural j=1;j<=5;j++) {
                for(Natural k=1;k<=5;k++) {
                    z[itok(k)] += D[ijktol(i,j,k,3,5)]*dx[itok(j)]*dy[itok(i)];
                }
            }
        }
    }
};

int main(int argc,char* argv[]){
    // Create a type shortcut
    using Optizelle::Rm;

    // Read in the name for the input file
    if(argc!=2) {
        std::cerr << "nw_sqp_exercise <parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string fname(argv[1]);

    // Generate an initial guess for the primal
    std::vector <double> x(5);
    x[itok(1)]=-1.8;
    x[itok(2)]=1.7;
    x[itok(3)]=1.9;
    x[itok(4)]=-0.8;
    x[itok(5)]=-0.8;

    // Generate an initial guess for the dual
    std::vector <double> y(3);

    // Create an optimization state
    Optizelle::EqualityConstrained <double,Rm,Rm>::State::t state(x,y);

    // Read the parameters from file
    Optizelle::json::EqualityConstrained <double,Optizelle::Rm,Optizelle::Rm>
        ::read(Optizelle::Messaging(),fname,state);
    
    // Create a bundle of functions
    Optizelle::EqualityConstrained <double,Rm,Rm>::Functions::t fns;
    fns.f.reset(new MyObj);
    fns.g.reset(new MyEq);
   
    // Solve the optimization problem
    Optizelle::EqualityConstrained <double,Rm,Rm>::Algorithms
        ::getMin(Optizelle::Messaging(),fns,state);

    // Print out the reason for convergence
    std::cout << "The algorithm converged due to: " <<
        Optizelle::StoppingCondition::to_string(state.opt_stop) << std::endl;

    // Print out the final answer
    std::cout << "The optimal point is:" << std::endl;
    for(Natural i=1;i<=5;i++) {
        if(i==1)
            std::cout << "[ ";
        else
            std::cout << "  ";
        std::cout << std::scientific << std::setprecision(16) << std::setw(23)
            << std::right << state.x[itok(i)];
        if(i==5)
            std::cout << " ]";
        else
            std::cout << " ;";
        std::cout << std::endl;
    }

    // Write out the final answer to file
    Optizelle::json::EqualityConstrained<double,Optizelle::Rm,Optizelle::Rm>
        ::write_restart(Optizelle::Messaging(),"solution.json",state);

    // Successful termination
    return EXIT_SUCCESS;
}
