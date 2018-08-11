// Run TCG on a system with negative curvature, but force it to exit on the
// first iteration

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.A = std::make_unique<Matrix>(Unit::Matrix <Real>::symm_nd(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.x_offset = std::make_unique <Vector> (
        Unit::Vector <Real>::basic(setup.m));

    // Trust-Region
    {
        // Make the TR radius large enough to not constrain ourselves until we
        // hit the negative curvature.
        setup.delta = 1e5;

        // Target solutions
        setup.iter_star = 1;
        setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;
        setup.x_star = std::make_unique <Vector> (*setup.b);
        auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
        X::scal((setup.delta-norm_b)/norm_b,*setup.x_star);

        // Tests
        Unit::reset_checks <Real,Rm> (setup);
        setup.check_iter = true;
        setup.check_stop = true;
        setup.check_tr = true;
        setup.check_cp = true;
        setup.check_sol = true;

        // Check the solver
        Unit::run_and_verify <Real,Rm> (setup);
    }

    // Line-search
    {
        // Blow out the trust-region
        setup.delta = std::numeric_limits<Real>::infinity();

        // Target solutions
        setup.iter_star = 1;
        setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;
        setup.x_star = std::make_unique <Vector> (*setup.b);

        // Tests
        Unit::reset_checks <Real,Rm> (setup);
        setup.check_iter = true;
        setup.check_stop = true;
        setup.check_cp = true;
        setup.check_sol = true;

        // Check the solver
        Unit::run_and_verify <Real,Rm> (setup);
    }

    // Declare success
    return EXIT_SUCCESS;
}
