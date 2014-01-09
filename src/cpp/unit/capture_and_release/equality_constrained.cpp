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

    // Create a different arbitrary vector in R^3
    std::vector <double> y(3);
    y[0] = 3.4; y[1] = 4.5; y[2]=5.6;

    // Create an unconstrained state based on this vector
    Optizelle::EqualityConstrained <double,Rm,Rm>::State::t state(x,y);

    // Do a release 
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart::X_Vectors xs;
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart::Y_Vectors ys;
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart::Reals reals;
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart::Naturals nats;
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart::Params params;
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart
        ::release(state,xs,ys,reals,nats,params);

    // Check that the state has empty slots for the variables
    CHECK(state.x.size()==0);
    CHECK(state.grad.size()==0);
    CHECK(state.dx.size()==0);
    CHECK(state.x_old.size()==0);
    CHECK(state.grad_old.size()==0);
    CHECK(state.dx_old.size()==0);
    CHECK(state.oldY.size()==0);
    CHECK(state.oldS.size()==0);
    CHECK(state.y.size()==0);
    CHECK(state.g_x.size()==0);
    CHECK(state.gpxdxn_p_gx.size()==0);
    CHECK(state.gpxdxt.size()==0);
    CHECK(state.dx_n.size()==0);
    CHECK(state.dx_ncp.size()==0);
    CHECK(state.dx_t.size()==0);
    CHECK(state.dx_t_uncorrected.size()==0);
    CHECK(state.dx_tcp_uncorrected.size()==0);
    CHECK(state.H_dxn.size()==0);
    CHECK(state.W_gradpHdxn.size()==0);
    CHECK(state.H_dxtuncorrected.size()==0);

    // Capture the state
    Optizelle::EqualityConstrained <double,Rm,Rm>::Restart
        ::capture(Optizelle::Messaging(),state,xs,ys,reals,nats,params);

    // Check that we actually have memory in these slots
    CHECK(state.x.size()>0);
    CHECK(state.grad.size()>0);
    CHECK(state.dx.size()>0);
    CHECK(state.x_old.size()>0);
    CHECK(state.grad_old.size()>0);
    CHECK(state.dx_old.size()>0);
    CHECK(state.y.size()>0);
    CHECK(state.g_x.size()>0);
    CHECK(state.gpxdxn_p_gx.size()>0);
    CHECK(state.gpxdxt.size()>0);
    CHECK(state.dx_n.size()>0);
    CHECK(state.dx_ncp.size()>0);
    CHECK(state.dx_t.size()>0);
    CHECK(state.dx_t_uncorrected.size()>0);
    CHECK(state.dx_tcp_uncorrected.size()>0);
    CHECK(state.H_dxn.size()>0);
    CHECK(state.W_gradpHdxn.size()>0);
    CHECK(state.H_dxtuncorrected.size()>0);

    // Check the relative error between the vector created above and the one
    // left in the state
    std::vector <double> residual(Rm <double>::init(y));
    Rm <double>::copy(y,residual);
    Rm <double>::axpy(-1.,state.y,residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(y,y)));
    CHECK(err < 1e-15);

    // Make sure we know we're successful
    return EXIT_SUCCESS;
}
