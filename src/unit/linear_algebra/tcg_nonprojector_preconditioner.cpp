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
        Unit::Matrix <Real>::symmetric(setup.m,25));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.proj_check = true;
    setup.iter_star = 1;
    setup.stop_star = Optizelle::TruncatedStop::NonProjector;
    setup.check_sol = false;
    setup.check_res = false;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
