// Run TCG with with a nonsymmetric preconditioner and show that we can detect
// it.  Note, this requires us to overorthogonalize. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::nonsymmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.orthog_storage_max = 3;

    // Target solutions
    setup.stop_star = Optizelle::TruncatedStop::NonSymmetricPreconditioner;

    // Tests
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
