// Run TCG on a system with negative curvature.  In one case, we have a
// trust-region with an offset.  In the other, we remove the trust-region like
// we would for a line-search method.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::mostly_dd_indef(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::elast(setup.m));
    setup.x_offset = std::make_unique <Vector> (
        Unit::Vector <Real>::basic(setup.m));

    // Trust-Region
    {
        // Make the TR radius large enough to not constrain ourselves until we
        // hit the negative curvature.
        setup.delta = 1e5;

        // Target solutions
        setup.iter_star = 5;
        setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;

        // Tests
        Unit::reset_checks <Real,Rm> (setup);
        setup.check_iter = true;
        setup.check_stop = true;
        setup.check_tr = true;

        // Check the solver
        Unit::run_and_verify <Real,Rm> (setup);
    }

    // Line-search
    {
        // Blow out the trust-region
        setup.delta = std::numeric_limits<Real>::infinity();

        // Target solutions
        setup.iter_star = 5;
        setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;

        // Tests
        Unit::reset_checks <Real,Rm> (setup);
        setup.check_iter = true;
        setup.check_stop = true;

        // Check the solver
        Unit::run_and_verify <Real,Rm> (setup);

        // Extract the solution that we generated.  This should be the
        // answer at the iteration prior to this.  Pull these in as the
        // answers and then rerun the example
        setup.iter_max = setup.iter-1;
        setup.x_star = std::make_unique <Vector> (*setup.x);
        setup.iter_star = setup.iter-1;
        setup.stop_star = Optizelle::TruncatedStop::MaxItersExceeded;

        // Tests
        Unit::reset_checks <Real,Rm> (setup);
        setup.check_stop = true;
        setup.check_sol = true;
        setup.check_iter = true;

        // Check the solver
        Unit::run_and_verify <Real,Rm> (setup);
    }

    // Declare success
    return EXIT_SUCCESS;
}
