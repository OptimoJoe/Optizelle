// This tests our ability to capture and release from the optimization state

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include "unit.h"
#include "restart.h"

int main() {

    // Create a messaging object
    Optizelle::Messaging msg;

    // Create some arbitrary vector in R^2
    std::vector <Real> x = {1.2,2.3};
    std::vector <Real> x0 = {2.3,1.2};
    
    // Create a different arbitrary vector in R^4
    std::vector <Real> z = {6.7,7.8,8.9,9.10};
    std::vector <Real> z0 = {9.10,8.9,7.8,6.7};

    // Create a state based on this vector
    //---State0---
    Optizelle::InequalityConstrained <Real,XX,ZZ>::State::t state(x,z);
    //---State1---

    // Read in some parameters
    std::string fname("blank.json");
    //---ReadJson0--- 
    Optizelle::json::InequalityConstrained <Real,XX,ZZ>::read(msg,fname,state);
    //---ReadJson1--- 
   
    // Create a bundle of functions
    //---Functions0---
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Functions::t fns;
    //---Functions1---
    fns.f.reset(new F);
    fns.h.reset(new H);

    // Do a null optimization
    //---Solver0---
    Optizelle::InequalityConstrained<Real,XX,ZZ>::Algorithms::getMin(
        msg,fns,state);
    //---Solver1---

    // Do a null optimization with a state manipulator 
    BlankManipulator <Optizelle::InequalityConstrained<Real,XX,ZZ> > smanip;
    //---SmanipSolver0---
    Optizelle::InequalityConstrained<Real,XX,ZZ>::Algorithms::getMin(
        msg,smanip,fns,state);
    //---SmanipSolver1---

    // Read and write the state to file
    fname = "restart.json";
    //---WriteReadRestart0---
    Optizelle::json::InequalityConstrained <Real,XX,ZZ>::write_restart(
        msg,fname,state);
    Optizelle::json::InequalityConstrained <Real,XX,ZZ>::read_restart(
        msg,fname,x,z,state);
    //---WriteReadRestart1---

    // Do a release 
    //---Release0---
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::X_Vectors xs;
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::Z_Vectors zs;
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::Reals reals;
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::Naturals nats;
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart::Params params;
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart
        ::release(state,xs,zs,reals,nats,params);
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
    CHECK(state.z.size()==0);
    CHECK(state.dz.size()==0);
    CHECK(state.h_x.size()==0);
    
    // Check that we have the correct number of vectors
    CHECK(xs.size() == 6);
    CHECK(zs.size() == 3);
    
    // Modify some vectors 
    xs.front().second = x0;
    zs.front().second = z0;

    // Capture the state
    //---Capture0---
    Optizelle::InequalityConstrained <Real,XX,ZZ>::Restart
        ::capture(msg,state,xs,zs,reals,nats,params);
    //---Capture1---

    // Check that we actually have memory in these slots
    CHECK(state.x.size()>0);
    CHECK(state.grad.size()>0);
    CHECK(state.dx.size()>0);
    CHECK(state.x_old.size()>0);
    CHECK(state.grad_old.size()>0);
    CHECK(state.dx_old.size()>0);
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
