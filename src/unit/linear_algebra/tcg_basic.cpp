// Run TCG with a preconditioner.  This verifies that TCG will solve a linear
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

    setup.x_star = std::make_unique <Vector> (Vector({
        9.75321699253918e-02,
        -2.05661763246187e-02,
        -6.58262181619375e-02,
        -3.55954161842943e-02,
        4.63421647030109e-03}));
    setup.iter_star = 5;
    setup.stop_star = Optizelle::TruncatedStop::RelativeErrorSmall;

    setup.check_sol = true;
    setup.check_iter = true;
    setup.check_res = true;
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
