// Test the quasinormal steps ability to not run when already feasible.  Here,
// we have two circles where the first has the center at (1,0)
// and the second has the center at (1,1).  The intersection occurs at
// (1+-sqrt(3)/2,1/2).

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess
    auto x = std::vector <Real> { 1.-std::sqrt(3.)/2., 0.5 };
    auto y = std::vector <Real> { 0., 0. };

    // Setup the test
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::Feasible;
    setup.dx_n_star.reset(new std::vector <Real> {0.0,0.0});
    setup.dx_ncp_star.reset(new std::vector <Real> {0.0,0.0});

    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success
    return EXIT_SUCCESS;
}
