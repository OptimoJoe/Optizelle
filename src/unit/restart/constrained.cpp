// This tests our ability to capture and release from the optimization state

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include "unit.h"
#include "restart.h"

int main() {

    // Create a messaging object
    auto msg = Optizelle::Messaging::stdout;

    // Create some arbitrary vector in R^2
    std::vector <Real> x = {1.2,2.3};
    std::vector <Real> x0 = {2.3,1.2};

    // Create a different arbitrary vector in R^3
    std::vector <Real> y = {3.4,4.5,5.6};
    std::vector <Real> y0 = {5.6,4.5,3.4};

    // Create a different arbitrary vector in R^4
    std::vector <Real> z = {6.7,7.8,8.9,9.10};
    std::vector <Real> z0 = {9.10,8.9,7.8,6.7};

    // Create a state based on this vector
    //---State0---
    Optizelle::Constrained <Real,XX,YY,ZZ>::State::t state(x,y,z);
    //---State1---

    // Read in some parameters
    std::string fname("blank.json");
    //---ReadJson0---
    Optizelle::json::Constrained <Real,XX,YY,ZZ>::read(fname,state);
    //---ReadJson1---

    // Create a bundle of functions
    //---Functions0---
    Optizelle::Constrained <Real,XX,YY,ZZ>::Functions::t fns;
    //---Functions1---
    fns.f.reset(new F);
    fns.g.reset(new G);
    fns.h.reset(new H);

    // Do a null optimization
    state.f_x = 1.0;
    //---Solver0---
    Optizelle::Constrained<Real,XX,YY,ZZ>::Algorithms::getMin(
        msg,fns,state);
    //---Solver1---

    // Do a null optimization with a state manipulator
    BlankManipulator <Optizelle::Constrained<Real,XX,YY,ZZ> > smanip;
    //---SmanipSolver0---
    Optizelle::Constrained<Real,XX,YY,ZZ>::Algorithms::getMin(
        msg,fns,state,smanip);
    //---SmanipSolver1---

    // Read and write the state to file
    fname = "restart.json";
    //---WriteReadRestart0---
    Optizelle::json::Constrained <Real,XX,YY,ZZ>::write_restart(
        fname,state);
    Optizelle::json::Constrained <Real,XX,YY,ZZ>::read_restart(
        fname,x,y,z,state);
    //---WriteReadRestart1---

    // Do a release
    //---Release0---
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::X_Vectors xs;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::Y_Vectors ys;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::Z_Vectors zs;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::Reals reals;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::Naturals nats;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart::Params params;
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart
        ::release(state,xs,ys,zs,reals,nats,params);
    //---Release1---

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
    CHECK(state.z.size()==0);
    CHECK(state.dz.size()==0);
    CHECK(state.h_x.size()==0);

    // Check that we have the correct number of vectors
    CHECK(xs.size() == 14);
    CHECK(ys.size() == 5);
    CHECK(zs.size() == 3);

    // Modify some vectors
    xs.front().second = x0;
    ys.front().second = y0;
    zs.front().second = z0;

    // Capture the state
    //---Capture0---
    Optizelle::Constrained <Real,XX,YY,ZZ>::Restart
        ::capture(state,xs,ys,zs,reals,nats,params);
    //---Capture1---

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
    CHECK(state.z.size()>0);
    CHECK(state.dz.size()>0);
    CHECK(state.h_x.size()>0);

    // Check the relative error between the vector created above and the one
    // left in the state.
    {
    std::vector <Real> residual(XX <Real>::init(x));
    XX <Real>::copy(x0,residual);
    XX <Real>::axpy(-1.,state.x,residual);
    Real err=std::sqrt(XX <Real>::innr(residual,residual))
        /(1.+sqrt(XX <Real>::innr(x0,x0)));
    CHECK(err < 1e-15);
    }

    {
    std::vector <Real> residual(YY <Real>::init(y));
    YY <Real>::copy(y0,residual);
    YY <Real>::axpy(-1.,state.y,residual);
    Real err=std::sqrt(YY <Real>::innr(residual,residual))
        /(1.+sqrt(YY <Real>::innr(y0,y0)));
    CHECK(err < 1e-15);
    }

    {
    std::vector <Real> residual(ZZ <Real>::init(z));
    ZZ <Real>::copy(z0,residual);
    ZZ <Real>::axpy(-1.,state.z,residual);
    Real err=std::sqrt(ZZ <Real>::innr(residual,residual))
        /(1.+sqrt(ZZ <Real>::innr(z0,z0)));
    CHECK(err < 1e-15);
    }

    // Make sure we know we're successful
    return EXIT_SUCCESS;
}
