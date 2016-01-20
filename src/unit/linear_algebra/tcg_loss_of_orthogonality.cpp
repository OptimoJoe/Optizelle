// Run TCG where we have a large enough poorly enough conditioned problem where
// we eventually lose orthogonality of the previous search directions. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.m=500;
    setup.eps=1e-13;
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b=std::make_unique <Vector> (Unit::Vector <Real>::alternate(setup.m));
    setup.orthog_storage_max=setup.m;
    setup.stop_star = Optizelle::TruncatedStop::LossOfOrthogonality;
    setup.check_iter=false;
    setup.check_sol=false;
    setup.check_res=false;
    setup.check_stop=true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
