// Run TCG with a offset that violates the safeguard 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
    setup.x_offset=std::make_unique<Vector>(Unit::Vector<Real>::ones(setup.m));
    X::scal(Real(-2.),*setup.x_offset);
    auto lb = Unit::Vector<Real>::zero(setup.m);
    auto x = Unit::Vector<Real>::ones(setup.m);
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,XX>>(
        Unit::Safeguard <Real,XX>::lower(x,lb));

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
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
