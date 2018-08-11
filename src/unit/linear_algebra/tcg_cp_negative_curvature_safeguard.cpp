// Run TCG on a system with negative curvature, but force it to exit on the
// first iteration.  In addition, put in a safeguard that cuts off the
// Cauchy-Point.  For reference, the Caucy Point is:
//
// 45585.961559216346
// -20585.918591763053
// -67831.200126303112
// -52712.789084326592
// 10869.517144297988
//
// which I found by just running the code.  As such, we just need to cut-off
// the second to last element.

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
    auto lb = Real(-10.);
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {0,0,0,1,0};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));
    setup.delta = 1e5;

    // Target solutions
    setup.iter_star = 1;
    setup.stop_star = Optizelle::TruncatedStop::NegativeCurvature;
    setup.x_star = std::make_unique <Vector> (*setup.b);
    auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
    X::scal((setup.delta-norm_b)/norm_b,*setup.x_star);

    // Tests
    setup.check_iter = true;
    setup.check_stop = true;
    setup.check_cp = true;
    setup.check_safeguard_alpha = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
