#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "unit.h"

// This tests our ability to capture and release from the optimization state

int main() {
    // Create a type shortcut
    using Optizelle::Rm;

    // Create some arbitrary vector in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    
    // Create a different arbitrary vector in R^4
    std::vector <double> z(4);
    z[0] = 6.7; z[1] = 7.8; z[2]=8.9; z[3]=9.10;

    // Create an unconstrained state based on this vector
    Optizelle::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Do a release 
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart::X_Vectors xs;
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart::Z_Vectors zs;
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart::Reals reals;
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart::Naturals nats;
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart::Params params;
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart
        ::release(state,xs,zs,reals,nats,params);

    // Check that the state has empty slots for the variables
    CHECK(state.x.size()==0);
    CHECK(state.grad.size()==0);
    CHECK(state.dx.size()==0);
    CHECK(state.x_old.size()==0);
    CHECK(state.grad_old.size()==0);
    CHECK(state.dx_old.size()==0);
    CHECK(state.oldY.size()==0);
    CHECK(state.oldS.size()==0);
    CHECK(state.z.size()==0);
    CHECK(state.h_x.size()==0);

    // Capture the state
    Optizelle::InequalityConstrained <double,Rm,Rm>::Restart
        ::capture(Optizelle::Messaging(),state,xs,zs,reals,nats,params);

    // Check that we actually have memory in these slots
    CHECK(state.x.size()>0);
    CHECK(state.grad.size()>0);
    CHECK(state.dx.size()>0);
    CHECK(state.x_old.size()>0);
    CHECK(state.grad_old.size()>0);
    CHECK(state.dx_old.size()>0);
    CHECK(state.z.size()>0);
    CHECK(state.h_x.size()>0);

    // Check the relative error between the vector created above and the one
    // left in the state.  We check x instead of z since z is actually set
    // internally and we don't use the variable provided by the user.
    std::vector <double> residual(Rm <double>::init(x));
    Rm <double>::copy(x,residual);
    Rm <double>::axpy(-1.,state.x,residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x,x)));
    CHECK(err < 1e-15);

    // Make sure we know we're successful
    return EXIT_SUCCESS;
}
