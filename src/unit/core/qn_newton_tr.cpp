// Test the quasinormal step where we cut off the Newton step and force a
// dogleg step on the intersection of two circles constraint.  At the same time,
// concoct things, so that the truncated Newton step works better than the
// dogleg step

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess
    auto x = std::vector <Real> { 0.9, 0.75 };
    auto y = std::vector <Real> { 0., 0. };

    // Setup the test
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));
    setup.delta = 0.8;

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::NewtonTrustRegion;

    // Set what tests we want
    setup.check_stop = true;
    setup.check_tr = true;
    setup.check_augsys = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success
    return EXIT_SUCCESS;
}
