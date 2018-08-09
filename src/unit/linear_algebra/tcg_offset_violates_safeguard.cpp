// Run TCG with a offset that violates the safeguard.  Specifically,
//
// x = (1,1,1,1,1)
// x_offset = (-2,-2,-2,-2,-2)
// w = (0,0,0,0,1)
// lb = 0
//
// Hence, <x,w> = 1 >= 0, but <x+x_offset,w> = -1 ~>= 0, which violates the
// bound.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
    setup.x_offset=std::make_unique<Vector>(Unit::Vector<Real>::ones(setup.m));
    X::scal(Real(-2.),*setup.x_offset);
    auto lb = Real(0.);
    auto x = Unit::Vector<Real>::ones(setup.m);
    auto w = std::vector <Real> {0,0,0,0,1};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));

    // Target solutions
    setup.x_star =std::make_unique<Vector>(Unit::Vector <Real>::basic(setup.m));
    X::zero(*setup.x_star);
    setup.iter_star = 0;
    setup.stop_star = Optizelle::TruncatedStop::OffsetViolatesSafeguard;

    // Tests
    setup.check_stop = true;
    setup.check_iter = true;
    setup.check_sol = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
