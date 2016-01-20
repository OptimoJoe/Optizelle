// Run TCG, but require too much accuracy.  Eventually, this should get a
// numerical instability that causes the objective metric to rise.  We need the
// objective metric, which is just the CG objective function, to always go down
// in order to guarantee model decrease for trust-region methods. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,XX> ();

    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.B = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,25));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.eps = 1e-24;

    setup.stop_star = Optizelle::TruncatedStop::ObjectiveIncrease;

    setup.check_stop = true;

    // Check the solver 
    Unit::run_and_verify <Real,XX> (setup);

    // Extract the solution that we generated.  This should be the answer at
    // the iteration prior to this.  Pull these in as the answers and then
    // rerun the example
    setup.iter_max = setup.iter-1;

    setup.x_star = std::make_unique <Vector> (*setup.x);
    setup.iter_star = setup.iter-1;
    setup.stop_star = Optizelle::TruncatedStop::MaxItersExceeded;

    setup.check_stop = true;
    setup.check_sol = true;
    setup.check_iter = true;
    setup.check_res = true;
    Unit::run_and_verify <Real,XX> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
