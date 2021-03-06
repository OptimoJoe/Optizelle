// Run TCG with a safeguard where the point prior to negative curvature is
// not safe and nor is the final point.  This means that we must retreat to the
// last safe point and do a safeguard step from there.  Specifically,
//
// x_offset = (0,0,0,0,0)
// w = (0,0,0,0,-1)
// lb = -31/40
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
// -1/2, -2/3, -3/4, -4/5, -13000
//
// As such, we require that <x + x_offset,w> >= -31/40, which means that we
// violate the safeguard on the last two steps.  This forces us to retreat
// to step 3 where <x + x_offset,w> is -3/4 and safeguard from there.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::mostly_dd_indef(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::elast(setup.m));
    auto lb = Real(-31./40.);
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {0,0,0,0,-1};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));
    setup.delta = 1e5;

    // Target solutions
    setup.iter_star = 5;
    setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;
    setup.safeguard_failed_star = 1;

    // Tests
    Unit::reset_checks <Real,Rm> (setup);
    setup.check_iter = true;
    setup.check_stop = true;
    setup.check_safeguard_failed = true;
    setup.check_safeguard_alpha = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
