// Run TCG where we increase the number of orthogonalization iterations.  This
// should allow us to solve the system whereas we have an error without this
// option turned on.

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.m=500;
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b=std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.eps=1e-13;
    setup.orthog_storage_max=setup.m;
    setup.orthog_iter_max=2;

    setup.stop_star = Optizelle::TruncatedStop::RelativeErrorSmall;

    setup.check_res = true;
    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
