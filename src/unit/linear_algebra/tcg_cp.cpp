// Run TCG, but limit ourselves to a single iteration.  This verifies that this
// solution is the Cauchy point.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.iter_max = 1;

    // Target solutions
    setup.iter_star = 1;
    setup.stop_star = Optizelle::TruncatedStop::MaxItersExceeded;

    // Tests
    setup.check_cp = true;
    setup.check_iter = true;
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
