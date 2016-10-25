// Run TCG with a safeguard where the point prior to negative curvature is
// safe, but the final point is not safe.  This forces a safeguard between the
// last point and the final point.  Specifically,
//
// x_offset = (0,0,0,0,0)
// w = (0,0,0,0,-1)
// lb = -10 
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
// As such, we require that <x + x_offset,w> >= -10, which means that we
// violate the safeguard on the final step.  However, the step prior to that is
// feasible, so we just do a linesearch as far as possible after detecting
// negative curvature.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::mostly_dd_indef(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::elast(setup.m));
    auto lb = Real(-10.); 
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {0,0,0,0,-1};
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
    setup.check_safeguard_failed = true;
    setup.check_safeguard_alpha = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
