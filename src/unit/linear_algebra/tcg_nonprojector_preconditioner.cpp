// Run TCG with with a nonprojector preconditioner and show that we can detect
// it. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,25));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.orthog_storage_max = 3;
    setup.check_B_projector = true;

    // Target solutions
    setup.stop_star = Optizelle::TruncatedStop::NonProjector;

    // Tests
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
