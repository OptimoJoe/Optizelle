// Run TCG with a small trust-region radius.  This verifies that the algorithm
// will not violate the trust-region radius and that the final solution will
// touch the radius.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.delta = Real(1e-3);
    setup.iter_max = 1;

    // Target solutions
    setup.iter_star = 1;
    setup.stop_star = Optizelle::TruncatedStop::TrustRegionViolated;

    // Tests
    setup.check_cp = true;
    setup.check_tr = true;
    setup.check_iter = true;
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
