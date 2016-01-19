// Run TCG with no preconditioners.  This verifies that TCG will solve a linear
// system.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.iter_star = 3;
    setup.stop_star = Optizelle::TruncatedStop::NonSymmetric;
    setup.check_sol = false;
    setup.check_res = false;
    setup.orthog_storage_max = 3;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
