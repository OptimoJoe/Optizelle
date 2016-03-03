// Test the quasinormal step where we cut off the Newton step and force a
// dogleg step on the intersection of two circles constraint.  At the same time,
// concoct things, so that the truncated Newton step works better than the
// dogleg step.  Finally, also cut off the Cauchy point with an inequality
// constraint.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "qn.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess 
    auto x = std::vector <Real> { 0.9, 0.75 };
    auto y = std::vector <Real> { 0., 0. };
  
    // Run an unrestricted test to get the Cauchy point
    auto cp = X::init(x);
    {
        // Setup the test 
        auto setup = Unit <Real>::QN(x,y);
        setup.g.reset(
            new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));

        // Set what tests we want
        setup.check_augsys = true;

        // Run the test
        Unit <Real>::run_and_verify(setup);

        // Grab the Cauchy point
        X::copy(setup.cp,cp);
    }

    // Now, cut off the Cauchy point
    {

        // Setup the test 
        auto setup = Unit <Real>::QN(x,y);
        setup.g.reset(
            new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));
        setup.h.reset(new Unit <Real>::Constraint::Box(
            {-10,-10.},{10.,0.8}));
        setup.delta = 0.8;

        // Set the targets
        setup.qn_stop_star = Optizelle::QuasinormalStop::NewtonTrustRegion;
        
        // Set what tests we want
        setup.check_stop = true;
        setup.check_tr = true;
        setup.check_augsys = true;

        // Run the test
        Unit <Real>::run_and_verify(setup);

        // Verify that we cut off the Cauchy point
        auto norm_r = Real(0.);
        auto norm_cp = Real(0.);
        std::tie(norm_r,norm_cp)=Unit <Real>::error(setup.cp,cp);
        CHECK(norm_r > 0.1);
    }

    // Declare success 
    return EXIT_SUCCESS;
}
