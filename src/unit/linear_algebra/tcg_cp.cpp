// Run TCG, but limit ourselves to a single iteration.  This verifies that this
// solution is the Cauchy point.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.iter_max = 1;
    setup.iter_star = 1;
    setup.check_cp = true;
    setup.check_sol = false;
    setup.check_res = false;
    setup.stop_star = Optizelle::TruncatedStop::MaxItersExceeded;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
