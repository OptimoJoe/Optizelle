// Run GMRES with the restart mechanism turned off.  This verifies that we can
// still converge even if we throw out Krylov vectors.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::gmres <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.x_star = std::make_unique <Vector> (Vector({
        6.71115708873876e-01,
        1.06789410922663e+00,
        -1.31466092004730e+00,
        5.25893325732259e-02,
        9.35912328721990e-01}));
    setup.rst_freq = 3;
    // Check our number of iterations.  Really, there's a question if this is
    // the right right number.  Basically, I just ran the code to create a
    // baseline and we want to know if this baseline changes.
    setup.iter_star = 227;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
