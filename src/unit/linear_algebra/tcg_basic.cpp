#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.x_star = std::make_unique <Vector> (Vector({
        9.56220622356691e-02,
        3.82201576243788e-03,
        -1.66128455554842e-01,
        -1.76933207992568e-01,
        1.79489652213403e-02}));
    setup.iter_star = 5;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
