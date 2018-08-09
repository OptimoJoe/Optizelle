// Forces the quasinormal step to exit on the first step with the Cauchy point

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess
    auto x = std::vector <Real> { 0., 0.};
    auto y = std::vector <Real> { 0. };

    // Setup the test
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::Linear);

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::CauchySolved;
    setup.dx_n_star.reset(new std::vector <Real> {0.5,0.5});
    setup.dx_ncp_star = std::make_unique <X_Vector> (*setup.dx_n_star);

    // Set what tests we want
    setup.check_stop = true;
    setup.check_feas = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success
    return EXIT_SUCCESS;
}
