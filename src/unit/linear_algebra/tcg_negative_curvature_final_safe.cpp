// Run TCG with a safeguard where the final point after detecting negative
// curvature is feasible, but the last point prior to that is not.  This means
// that we can take the final point without truncating due to the safeguard.
// Specifically,
//
// x_offset = (0,0,0,0,0)
// w = (4,4,4,1,1)
// lb = -0.1
//
// The sequence of solutions produced by CG is approximately
//
// 0 0 0 0 0
// 0 0 0 0 1/2
// 0 0 0 -1/3 2/3
// 0 0 1/4 -1/2 3/4
// 0 -1/5 2/5 -3/5 4/5
// 67000 -53000 40000 -27000 13000
//
// Hence, the sequence <x + x_offset,w> is
//
// 0, 1/2, 1/3, 1/2, -2/5, 82000
//
// As such, we require that <x + x_offset,w> >= -1/5, which means that we
// violate the safeguard on the second to last step, but satisfy it when we
// take a step after detecting negative curvature.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::mostly_dd_indef(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::elast(setup.m));
    auto lb = Real(-0.1);
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {4,4,1,1,1};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));
    setup.delta = 1e5;

    // Target solutions
    setup.iter_star = 5;
    setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;
    setup.safeguard_failed_star = 0;

    // Tests
    Unit::reset_checks <Real,Rm> (setup);
    setup.check_iter = true;
    setup.check_stop = true;
    setup.check_tr = true;
    setup.check_safeguard_failed = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
