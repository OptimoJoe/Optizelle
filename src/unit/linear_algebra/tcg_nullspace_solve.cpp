// Run TCG with a projection operator.  This verifies that with a projection
// operator, we'll converge in a smaller number of iterations and that
// projections that are singular don't break anything too badly.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::project_2(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::sum_2(setup.m,0));

    // Target solutions
    setup.x_star = std::make_unique <Vector> (Vector({
        1.0,
        1.0,
        0.0,
        0.0,
        0.0}));
    setup.iter_star = 2;
    setup.stop_star = Optizelle::TruncatedStop::RelativeErrorSmall;

    // Tests
    setup.check_sol = true;
    setup.check_iter = true;
    setup.check_res = true;
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
