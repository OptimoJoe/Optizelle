// Test the quasinormal steps ability to exit on the Newton step when the
// Cauchy point violates the inequality constraint. 

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess 
    auto x = std::vector <Real> { 0., 0. };
    auto y = std::vector <Real> { 0., 0. };

    // Setup the test 
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));

    // Create an inequality constraint to cut off the Cauchy point 
    auto c = setup.zeta * setup.gamma;
    setup.h.reset(new Unit <Real>::Constraint::Box(
        {-1.,-1.},{0.1/c,10.}));

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::Newton;
    setup.dx_n_star.reset(new std::vector <Real> {0.0,0.5});
    setup.dx_ncp_star.reset(new std::vector <Real> {0.1,0.1});
    X::axpy(Real(-1.),x,*setup.dx_n_star);
    
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
