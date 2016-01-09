#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // On the first iteration, CG will move in the steepest descent direction.
    // If we place the offset in the same direction, we know that it'll still
    // move in that direction, but less so.  Here, we want 40% of the motion
    // to come from the offset and 60% from the steepest descent direction.
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.delta = Real(0.1);
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    auto norm_b = std::sqrt(X::innr(*setup.b,*setup.b));
    setup.x_offset=std::make_unique<Vector>(Unit::Vector<Real>::basic(setup.m));
    X::scal(Real(0.4)/norm_b*setup.delta,*setup.x_offset);
    setup.x_star =std::make_unique<Vector>(Unit::Vector <Real>::basic(setup.m));
    X::scal(Real(0.6)/norm_b*setup.delta,*setup.x_star);
    setup.iter_max = 1;
    setup.iter_star = 1;
    setup.check_tr = true;
    setup.check_res = false;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
