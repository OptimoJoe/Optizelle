// Run TCG with a small trust-region radius and an invalid offset.  

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // Problem setup
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.delta = Real(0.1);
    auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
    setup.x_offset=std::make_unique<Vector>(Unit::Vector<Real>::basic(setup.m));
    X::scal(Real(2.0)/norm_b*setup.delta,*setup.x_offset);

    // Target solutions
    setup.x_star =std::make_unique<Vector>(Unit::Vector <Real>::basic(setup.m));
    X::zero(*setup.x_star);
    setup.iter_star = 0;
    setup.stop_star = Optizelle::TruncatedStop::InvalidTrustRegionOffset;

    // Tests
    setup.check_stop = true;
    setup.check_iter = true;
    setup.check_sol = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
