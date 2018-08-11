// Run TCG where we have a large enough poorly enough conditioned problem where
// we eventually lose orthogonality of the previous search directions.  Then,
// increase the number of orthogonalization iterations and show we can solve
// the problem.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup
    setup.m=500;
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b=std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.orthog_storage_max=setup.m;
    setup.eps=1e-13;

    // Target solutions
    setup.stop_star = Optizelle::TruncatedStop::LossOfOrthogonality;

    // Tests
    setup.check_stop = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Change the problem setup to regain orthogonality
    setup.orthog_iter_max=2;

    // Target solutions
    setup.stop_star = Optizelle::TruncatedStop::RelativeErrorSmall;

    // Tests
    setup.check_res = true;
    setup.check_stop = true;

    // Check the solver
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
