// This tests our ability to capture and release from the optimization state

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
#include "unit.h"

// Create some type shortcuts
typedef double Real;
template <typename Real> using XX = Optizelle::Rm <Real>;
template <typename Real> using YY = Optizelle::Rm <Real>;
template <typename Real> using ZZ = Optizelle::Rm <Real>;

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
        
        template <typename Real>
        struct Serialization <Real,YY> {
            static std::string serialize(
                typename YY <Real>::Vector const & x
            ) {
                return Optizelle::json::Serialization
                    <Real,Optizelle::Rm>::serialize(x);
            }
            static typename YY <Real>::Vector deserialize(
                typename YY <Real>::Vector const & x,
                std::string const & x_json
            ) {
                return std::move(Optizelle::json::Serialization
                    <Real,Optizelle::Rm>::deserialize(x,x_json));
            }
        };
        
        template <typename Real>
        struct Serialization <Real,ZZ> {
            static std::string serialize(
                typename ZZ <Real>::Vector const & x
            ) {
                return Optizelle::json::Serialization
                    <Real,Optizelle::Rm>::serialize(x);
            }
            static typename ZZ <Real>::Vector deserialize(
                typename ZZ <Real>::Vector const & x,
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

    // Create a different arbitrary vector in R^3
    std::vector <Real> y = {3.4,4.5,5.6};
    std::vector <Real> y0 = {5.6,4.5,3.4};
    
    // Create a different arbitrary vector in R^4
    std::vector <Real> z = {6.7,7.8,8.9,9.10};
    std::vector <Real> z0 = {9.10,8.9,7.8,6.7};

    // Create a state based on this vector
    Optizelle::Constrained <Real,XX,YY,ZZ>::State::t state(x,y,z);

    // Read and write the state to file
    std::string fname("restart.json");
    //---WriteReadRestart0---
    Optizelle::json::Constrained <Real,XX,YY,ZZ>::write_restart(
        msg,fname,state);
    Optizelle::json::Constrained <Real,XX,YY,ZZ>::read_restart(
        msg,fname,x,y,z,state);
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
        ::capture(msg,state,xs,ys,zs,reals,nats,params);
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
