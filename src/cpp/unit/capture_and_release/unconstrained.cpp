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

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <double,Rm>::State::t state(x);

    // Do a release 
    Optizelle::Unconstrained <double,Rm>::X_Vectors xs;
    Optizelle::Unconstrained <double,Rm>::Reals reals;
    Optizelle::Unconstrained <double,Rm>::Nats nats;
    Optizelle::Unconstrained <double,Rm>::Params params;
    Optizelle::Unconstrained <double,Rm>::Restart
        ::release(state,xs,reals,nats,params);

    // Check that the state has empty slots for the variables
    CHECK(state.x.size()==0);
    CHECK(state.grad.size()==0);
    CHECK(state.dx.size()==0);
    CHECK(state.x_old.size()==0);
    CHECK(state.grad_old.size()==0);
    CHECK(state.dx_old.size()==0);
    CHECK(state.oldY.size()==0);
    CHECK(state.oldS.size()==0);

    // Capture the state
    Optizelle::Unconstrained <double,Rm>::Restart
        ::capture(Optizelle::Messaging(),state,xs,reals,nats,params);

    // Check that we actually have memory in these slots
    CHECK(state.x.size()>0);
    CHECK(state.grad.size()>0);
    CHECK(state.dx.size()>0);
    CHECK(state.x_old.size()>0);
    CHECK(state.grad_old.size()>0);
    CHECK(state.dx_old.size()>0);

    // Check the relative error between the vector created above and the one
    // left in the state
    std::vector <double> residual(Rm <double>::init(x));
    Rm <double>::copy(x,residual);
    Rm <double>::axpy(-1.,state.x,residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x,x)));
    CHECK(err < 1e-15);

    // Make sure we know we're successful
    return EXIT_SUCCESS;
}
