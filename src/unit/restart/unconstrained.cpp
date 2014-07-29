// This tests our ability to capture and release from the optimization state

//---Import0---
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
//---Import1---
#include "unit.h"
#include "restart.h"

template <typename Real> using WW = Optizelle::Rm <Real>;

//---Serialization0---
namespace Optizelle {
    namespace json {
        template <>
        struct Serialization <Real,WW> {
            static std::string serialize(
                typename WW <Real>::Vector const & x,
                std::string const & name,
                Natural const & iter
            ) { throw; }
            static typename WW <Real>::Vector deserialize(
                typename WW <Real>::Vector const & x,
                std::string const & x_json
            ) { throw; }
        };
    }
}
//---Serialization1---

int main() {

    // Create a messaging object
    Optizelle::Messaging msg;
    
    // Create some arbitrary vector in R^2
    std::vector <Real> x = {1.2,2.3};
    std::vector <Real> x0 = {2.3,1.2};

    // Create an unconstrained state based on this vector
    //---State0---
    Optizelle::Unconstrained <Real,XX>::State::t state(x);
    //---State1---

    // Read in some parameters
    std::string fname("blank.json");
    //---ReadJson0--- 
    Optizelle::json::Unconstrained <Real,XX>::read(msg,fname,state);
    //---ReadJson1--- 
   
    // Create a bundle of functions
    //---Functions0---
    Optizelle::Unconstrained <Real,XX>::Functions::t fns;
    //---Functions1---
    fns.f.reset(new F);

    // Do a null optimization
    //---Solver0---
    Optizelle::Unconstrained<Real,XX>::Algorithms::getMin(
        msg,fns,state);
    //---Solver1---

    // Do a null optimization with a state manipulator 
    BlankManipulator <Optizelle::Unconstrained<Real,XX> > smanip;
    //---SmanipSolver0---
    Optizelle::Unconstrained<Real,XX>::Algorithms::getMin(
        msg,fns,state,smanip);
    //---SmanipSolver1---

    // Read and write the state to file
    fname = "restart.json";
    //---WriteReadRestart0---
    Optizelle::json::Unconstrained <Real,XX>::write_restart(
        msg,fname,state);
    Optizelle::json::Unconstrained <Real,XX>::read_restart(
        msg,fname,x,state);
    //---WriteReadRestart1---

    // Do a release 
    //---Release0---
    Optizelle::Unconstrained <Real,XX>::Restart::X_Vectors xs;
    Optizelle::Unconstrained <Real,XX>::Restart::Reals reals;
    Optizelle::Unconstrained <Real,XX>::Restart::Naturals nats;
    Optizelle::Unconstrained <Real,XX>::Restart::Params params;
    Optizelle::Unconstrained <Real,XX>::Restart
        ::release(state,xs,reals,nats,params);
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

    // Check that we have the correct number of vectors
    CHECK(xs.size() == 6);

    // Modify some vectors 
    xs.front().second = x0;

    // Capture the state
    //---Capture0---
    Optizelle::Unconstrained <Real,XX>::Restart
        ::capture(msg,state,xs,reals,nats,params);
    //---Capture1---

    // Check that we actually have memory in these slots
    CHECK(state.x.size()>0);
    CHECK(state.grad.size()>0);
    CHECK(state.dx.size()>0);
    CHECK(state.x_old.size()>0);
    CHECK(state.grad_old.size()>0);
    CHECK(state.dx_old.size()>0);

    // Check the relative error between the vector created above and the one
    // left in the state
    {
    std::vector <Real> residual(XX <Real>::init(x));
    XX <Real>::copy(x0,residual);
    XX <Real>::axpy(-1.,state.x,residual);
    Real err=std::sqrt(XX <Real>::innr(residual,residual))
        /(1.+sqrt(XX <Real>::innr(x0,x0)));
    CHECK(err < 1e-15);
    }

    // Make sure we know we're successful
    return EXIT_SUCCESS;
}
