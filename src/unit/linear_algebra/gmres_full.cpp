// Run GMRES with both preconditioners set to something random.  We want to
// verify that we'll eventually find a solution.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::gmres <Real,Rm> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.B_left = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,30));
    setup.B_right = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,55));

    setup.x_star = std::make_unique <Vector> (Vector({
        6.71115708873876e-01,
        1.06789410922663e+00,
        -1.31466092004730e+00,
        5.25893325732259e-02,
        9.35912328721990e-01}));
    setup.iter_star = 5;

    setup.check_sol=true;
    setup.check_iter=true;
    setup.check_res=true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
