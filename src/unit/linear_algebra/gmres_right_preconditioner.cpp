// Run GMRES with a perfect right preconditioner.  Note, the operator here is
// singular, but this allows us to verify that we still converge in a single
// iteration.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::gmres <Real,Rm> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::diagonal(setup.m));
    setup.b = std::make_unique <Vector>(Unit::Vector<Real>::alternate(setup.m));
    setup.B_right = std::make_unique <Matrix>(
        Unit::Matrix <Real>::diagonal_inv(setup.m));

    setup.x_star = std::make_unique <Vector> (Vector({
        1.0,
        0.,
        0.5,
        0.,
        1./3.}));
    setup.iter_star = 1;

    setup.check_sol=true;
    setup.check_iter=true;
    setup.check_res=true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
