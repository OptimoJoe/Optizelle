#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::project_2(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::sum_2(setup.m,0));
    setup.x_star = std::make_unique <Vector> (Vector({
        1.0,
        1.0,
        0.0,
        0.0,
        0.0}));
    setup.iter_star = 2;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
