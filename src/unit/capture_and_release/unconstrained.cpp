// This tests our ability to capture and release from the optimization state

//---Import0---
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
//---Import1---
#include "unit.h"

// Create some type shortcuts
typedef double Real;
template <typename Real> using XX = Optizelle::Rm <Real>;

// There's some irritating bug in GCC that forces us to do this 
namespace Optizelle {
    namespace json {
        template <typename Real>
        struct Serialization <Real,XX> {
            static std::string serialize(
                typename XX <Real>::Vector const & x
            ) {
                return Optizelle::json::Serialization
                    <Real,Optizelle::Rm>::serialize(x);
            }
            static typename XX <Real>::Vector deserialize(
                typename XX <Real>::Vector const & x,
                std::string const & x_json
            ) {
                return std::move(Optizelle::json::Serialization
                    <Real,Optizelle::Rm>::deserialize(x,x_json));
            }
        };
    }
}

int main() {

    // Create a messaging object
    Optizelle::Messaging msg;
    
    // Create some arbitrary vector in R^2
    std::vector <Real> x = {1.2,2.3};
    std::vector <Real> x0 = {2.3,1.2};

    // Create an unconstrained state based on this vector
    Optizelle::Unconstrained <Real,XX>::State::t state(x);

    // Read and write the state to file
    std::string fname("restart.json");
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
