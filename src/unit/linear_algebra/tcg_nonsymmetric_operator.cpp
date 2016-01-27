// Run TCG with with a nonsymmetric operator and show that we can detect
// it.  Note, this requires us to overorthogonalize.  Note, this test is a
// little touchy since a nonsymmetric operator may cause the negative curvature
// or objective increase conditions to trigger first.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,4));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.orthog_storage_max = 3;

    // Target solutions
    setup.stop_star = Optizelle::TruncatedStop::NonSymmetricOperator;

    // Tests
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
