#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// This tests our ability to capture and release from the optimization state

BOOST_AUTO_TEST_SUITE(capture_and_release)

BOOST_AUTO_TEST_CASE(unconstrained) {
    // Create a type shortcut
    using peopt::Rm;

    // Create some arbitrary vector in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;

    // Create an unconstrained state based on this vector
    peopt::Unconstrained <double,Rm>::State::t state(x);

    // Do a release 
    peopt::Unconstrained <double,Rm>::X_Vectors xs;
    peopt::Unconstrained <double,Rm>::Reals reals;
    peopt::Unconstrained <double,Rm>::Nats nats;
    peopt::Unconstrained <double,Rm>::Params params;
    peopt::Unconstrained <double,Rm>::Restart
        ::release(state,xs,reals,nats,params);

    // Check that the state has empty slots for the variables
    BOOST_CHECK(state.x.size()==0);
    BOOST_CHECK(state.g.size()==0);
    BOOST_CHECK(state.dx.size()==0);
    BOOST_CHECK(state.x_old.size()==0);
    BOOST_CHECK(state.g_old.size()==0);
    BOOST_CHECK(state.dx_old.size()==0);
    BOOST_CHECK(state.oldY.size()==0);
    BOOST_CHECK(state.oldS.size()==0);

    // Capture the state
    peopt::Unconstrained <double,Rm>::Restart
        ::capture(peopt::Messaging(),state,xs,reals,nats,params);

    // Check that we actually have memory in these slots
    BOOST_CHECK(state.x.size()==1);
    BOOST_CHECK(state.g.size()==1);
    BOOST_CHECK(state.dx.size()==1);
    BOOST_CHECK(state.x_old.size()==1);
    BOOST_CHECK(state.g_old.size()==1);
    BOOST_CHECK(state.dx_old.size()==1);

    // Check the relative error between the vector created above and the one
    // left in the state
    std::vector <double> residual;
    Rm <double>::init(x,residual);
    Rm <double>::copy(x,residual);
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x,x)));
    BOOST_CHECK(err < 1e-15);
}

BOOST_AUTO_TEST_CASE(equality_constrained) {
    // Create a type shortcut
    using peopt::Rm;

    // Create some arbitrary vector in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;

    // Create a different arbitrary vector in R^3
    std::vector <double> y(3);
    y[0] = 3.4; y[1] = 4.5; y[2]=5.6;

    // Create an unconstrained state based on this vector
    peopt::EqualityConstrained <double,Rm,Rm>::State::t state(x,y);

    // Do a release 
    peopt::EqualityConstrained <double,Rm,Rm>::X_Vectors xs;
    peopt::EqualityConstrained <double,Rm,Rm>::Y_Vectors ys;
    peopt::EqualityConstrained <double,Rm,Rm>::Reals reals;
    peopt::EqualityConstrained <double,Rm,Rm>::Nats nats;
    peopt::EqualityConstrained <double,Rm,Rm>::Params params;
    peopt::EqualityConstrained <double,Rm,Rm>::Restart
        ::release(state,xs,ys,reals,nats,params);

    // Check that the state has empty slots for the variables
    BOOST_CHECK(state.x.size()==0);
    BOOST_CHECK(state.g.size()==0);
    BOOST_CHECK(state.dx.size()==0);
    BOOST_CHECK(state.x_old.size()==0);
    BOOST_CHECK(state.g_old.size()==0);
    BOOST_CHECK(state.dx_old.size()==0);
    BOOST_CHECK(state.oldY.size()==0);
    BOOST_CHECK(state.oldS.size()==0);
    BOOST_CHECK(state.y.size()==0);
    BOOST_CHECK(state.g_x.size()==0);
    BOOST_CHECK(state.gpxdxn_p_gx.size()==0);
    BOOST_CHECK(state.dx_n.size()==0);
    BOOST_CHECK(state.dx_ncp.size()==0);
    BOOST_CHECK(state.dx_t.size()==0);
    BOOST_CHECK(state.dx_t_uncorrected.size()==0);
    BOOST_CHECK(state.dx_tcp_uncorrected.size()==0);
    BOOST_CHECK(state.H_dxn.size()==0);
    BOOST_CHECK(state.W_gpHdxn.size()==0);
    BOOST_CHECK(state.H_dxtuncorrected.size()==0);

    // Capture the state
    peopt::EqualityConstrained <double,Rm,Rm>::Restart
        ::capture(peopt::Messaging(),state,xs,ys,reals,nats,params);

    // Check that we actually have memory in these slots
    BOOST_CHECK(state.x.size()==1);
    BOOST_CHECK(state.g.size()==1);
    BOOST_CHECK(state.dx.size()==1);
    BOOST_CHECK(state.x_old.size()==1);
    BOOST_CHECK(state.g_old.size()==1);
    BOOST_CHECK(state.dx_old.size()==1);
    BOOST_CHECK(state.y.size()==1);
    BOOST_CHECK(state.g_x.size()==1);
    BOOST_CHECK(state.gpxdxn_p_gx.size()==1);
    BOOST_CHECK(state.dx_n.size()==1);
    BOOST_CHECK(state.dx_ncp.size()==1);
    BOOST_CHECK(state.dx_t.size()==1);
    BOOST_CHECK(state.dx_t_uncorrected.size()==1);
    BOOST_CHECK(state.dx_tcp_uncorrected.size()==1);
    BOOST_CHECK(state.H_dxn.size()==1);
    BOOST_CHECK(state.W_gpHdxn.size()==1);
    BOOST_CHECK(state.H_dxtuncorrected.size()==1);

    // Check the relative error between the vector created above and the one
    // left in the state
    std::vector <double> residual;
    Rm <double>::init(y,residual);
    Rm <double>::copy(y,residual);
    Rm <double>::axpy(-1,state.y.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(y,y)));
    BOOST_CHECK(err < 1e-15);
}

BOOST_AUTO_TEST_CASE(inequality_constrained) {
    // Create a type shortcut
    using peopt::Rm;

    // Create some arbitrary vector in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    
    // Create a different arbitrary vector in R^4
    std::vector <double> z(4);
    z[0] = 6.7; z[1] = 7.8; z[2]=8.9; z[3]=9.10;

    // Create an unconstrained state based on this vector
    peopt::InequalityConstrained <double,Rm,Rm>::State::t state(x,z);

    // Do a release 
    peopt::InequalityConstrained <double,Rm,Rm>::X_Vectors xs;
    peopt::InequalityConstrained <double,Rm,Rm>::Z_Vectors zs;
    peopt::InequalityConstrained <double,Rm,Rm>::Reals reals;
    peopt::InequalityConstrained <double,Rm,Rm>::Nats nats;
    peopt::InequalityConstrained <double,Rm,Rm>::Params params;
    peopt::InequalityConstrained <double,Rm,Rm>::Restart
        ::release(state,xs,zs,reals,nats,params);

    // Check that the state has empty slots for the variables
    BOOST_CHECK(state.x.size()==0);
    BOOST_CHECK(state.g.size()==0);
    BOOST_CHECK(state.dx.size()==0);
    BOOST_CHECK(state.x_old.size()==0);
    BOOST_CHECK(state.g_old.size()==0);
    BOOST_CHECK(state.dx_old.size()==0);
    BOOST_CHECK(state.oldY.size()==0);
    BOOST_CHECK(state.oldS.size()==0);
    BOOST_CHECK(state.z.size()==0);
    BOOST_CHECK(state.h_x.size()==0);
    BOOST_CHECK(state.g_orig.size()==0);
    BOOST_CHECK(state.g_schur.size()==0);
    BOOST_CHECK(state.g_lag.size()==0);

    // Capture the state
    peopt::InequalityConstrained <double,Rm,Rm>::Restart
        ::capture(peopt::Messaging(),state,xs,zs,reals,nats,params);

    // Check that we actually have memory in these slots
    BOOST_CHECK(state.x.size()==1);
    BOOST_CHECK(state.g.size()==1);
    BOOST_CHECK(state.dx.size()==1);
    BOOST_CHECK(state.x_old.size()==1);
    BOOST_CHECK(state.g_old.size()==1);
    BOOST_CHECK(state.dx_old.size()==1);
    BOOST_CHECK(state.z.size()==1);
    BOOST_CHECK(state.h_x.size()==1);
    BOOST_CHECK(state.g_orig.size()==1);
    BOOST_CHECK(state.g_schur.size()==1);
    BOOST_CHECK(state.g_lag.size()==1);

    // Check the relative error between the vector created above and the one
    // left in the state.  We check x instead of z since z is actually set
    // internally and we don't use the variable provided by the user.
    std::vector <double> residual;
    Rm <double>::init(x,residual);
    Rm <double>::copy(x,residual);
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x,x)));
    BOOST_CHECK(err < 1e-15);
}

BOOST_AUTO_TEST_CASE(constrained) {
    // Create a type shortcut
    using peopt::Rm;

    // Create some arbitrary vector in R^2
    std::vector <double> x(2);
    x[0] = 1.2; x[1] = 2.3;
    
    // Create a different arbitrary vector in R^3
    std::vector <double> y(3);
    y[0] = 3.4; y[1] = 4.5; y[2]=5.6;

    // Create a different arbitrary vector in R^4
    std::vector <double> z(4);
    z[0] = 6.7; z[1] = 7.8; z[2]=8.9; z[3]=9.10;

    // Create an unconstrained state based on this vector
    peopt::Constrained <double,Rm,Rm,Rm>::State::t state(x,y,z);

    // Do a release 
    peopt::Constrained <double,Rm,Rm,Rm>::X_Vectors xs;
    peopt::Constrained <double,Rm,Rm,Rm>::Y_Vectors ys;
    peopt::Constrained <double,Rm,Rm,Rm>::Z_Vectors zs;
    peopt::Constrained <double,Rm,Rm,Rm>::Reals reals;
    peopt::Constrained <double,Rm,Rm,Rm>::Nats nats;
    peopt::Constrained <double,Rm,Rm,Rm>::Params params;
    peopt::Constrained <double,Rm,Rm,Rm>::Restart
        ::release(state,xs,ys,zs,reals,nats,params);

    // Check that the state has empty slots for the variables
    BOOST_CHECK(state.x.size()==0);
    BOOST_CHECK(state.g.size()==0);
    BOOST_CHECK(state.dx.size()==0);
    BOOST_CHECK(state.x_old.size()==0);
    BOOST_CHECK(state.g_old.size()==0);
    BOOST_CHECK(state.dx_old.size()==0);
    BOOST_CHECK(state.oldY.size()==0);
    BOOST_CHECK(state.oldS.size()==0);
    BOOST_CHECK(state.y.size()==0);
    BOOST_CHECK(state.g_x.size()==0);
    BOOST_CHECK(state.gpxdxn_p_gx.size()==0);
    BOOST_CHECK(state.dx_n.size()==0);
    BOOST_CHECK(state.dx_ncp.size()==0);
    BOOST_CHECK(state.dx_t.size()==0);
    BOOST_CHECK(state.dx_t_uncorrected.size()==0);
    BOOST_CHECK(state.dx_tcp_uncorrected.size()==0);
    BOOST_CHECK(state.H_dxn.size()==0);
    BOOST_CHECK(state.W_gpHdxn.size()==0);
    BOOST_CHECK(state.H_dxtuncorrected.size()==0);
    BOOST_CHECK(state.z.size()==0);
    BOOST_CHECK(state.h_x.size()==0);
    BOOST_CHECK(state.g_orig.size()==0);
    BOOST_CHECK(state.g_schur.size()==0);
    BOOST_CHECK(state.g_lag.size()==0);

    // Capture the state
    peopt::Constrained <double,Rm,Rm,Rm>::Restart
        ::capture(peopt::Messaging(),state,xs,ys,zs,reals,nats,params);

    // Check that we actually have memory in these slots
    BOOST_CHECK(state.x.size()==1);
    BOOST_CHECK(state.g.size()==1);
    BOOST_CHECK(state.dx.size()==1);
    BOOST_CHECK(state.x_old.size()==1);
    BOOST_CHECK(state.g_old.size()==1);
    BOOST_CHECK(state.dx_old.size()==1);
    BOOST_CHECK(state.y.size()==1);
    BOOST_CHECK(state.g_x.size()==1);
    BOOST_CHECK(state.gpxdxn_p_gx.size()==1);
    BOOST_CHECK(state.dx_n.size()==1);
    BOOST_CHECK(state.dx_ncp.size()==1);
    BOOST_CHECK(state.dx_t.size()==1);
    BOOST_CHECK(state.dx_t_uncorrected.size()==1);
    BOOST_CHECK(state.dx_tcp_uncorrected.size()==1);
    BOOST_CHECK(state.H_dxn.size()==1);
    BOOST_CHECK(state.W_gpHdxn.size()==1);
    BOOST_CHECK(state.H_dxtuncorrected.size()==1);
    BOOST_CHECK(state.z.size()==1);
    BOOST_CHECK(state.h_x.size()==1);
    BOOST_CHECK(state.g_orig.size()==1);
    BOOST_CHECK(state.g_schur.size()==1);
    BOOST_CHECK(state.g_lag.size()==1);

    // Check the relative error between the vector created above and the one
    // left in the state.  
    std::vector <double> residual;
    Rm <double>::init(x,residual);
    Rm <double>::copy(x,residual);
    Rm <double>::axpy(-1,state.x.front(),residual);
    double err=std::sqrt(Rm <double>::innr(residual,residual))
        /(1+sqrt(Rm <double>::innr(x,x)));
    BOOST_CHECK(err < 1e-15);
}

BOOST_AUTO_TEST_SUITE_END()
