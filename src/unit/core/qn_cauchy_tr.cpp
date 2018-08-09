// Verifies that we can cut off the Cauchy point with the trust-region.  In this
// setup, we know the actual Cauchy-point would be (0.5,0.5), so we just
// cut it halfway to be (0.25,0.25).

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
    // || (0.25,0.25) || = sqrt(1/8)
    setup.delta = std::sqrt(Real(1./8.))/setup.zeta;

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::CauchyTrustRegion;
    setup.dx_n_star.reset(new std::vector <Real> {0.25,0.25});
    setup.dx_ncp_star = std::make_unique <X_Vector> (*setup.dx_n_star);

    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;
    setup.check_tr = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success
    return EXIT_SUCCESS;
}
