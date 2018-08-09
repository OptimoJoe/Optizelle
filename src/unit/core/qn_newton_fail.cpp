// Test the quasinormal step where the derivative of our equality constraint is
// rank-deficient, which causes our augmented system solve for the Newton step
// to fail.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess
    auto x = std::vector <Real> { 1., 0.75 };
    auto y = std::vector <Real> { 0., 0. };

    // Setup the test
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::NewtonFailed;
    setup.dx_ncp_star.reset(new std::vector <Real> {0,0.075});
    setup.dx_n_star.reset(new std::vector <Real> (*setup.dx_ncp_star));

    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;
    setup.check_augsys = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success
    return EXIT_SUCCESS;
}
