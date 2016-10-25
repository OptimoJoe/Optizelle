// Run TCG, but limit ourselves to a single iteration, which forces to exit
// after the grabbing the Cauchy point, which happens to be
//
// 0.043401461872297958
// -0.019599432156455165
// -0.064580698647975734
// -0.050186768631698325
// 0.010348645016420882 
//
// which can be found by running the problem.  In addition, though, we add a
// safeguard to cut off the second to last element, which should mean that the
// Cauchy point is cut off. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::symmetric(setup.m,0));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::basic(setup.m));
    setup.iter_max = 1;
    auto lb = Real(-0.01);
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {0,0,0,1,0};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));

    // Target solutions
    setup.iter_star = 1;
    setup.stop_star = Optizelle::TruncatedStop::MaxItersExceeded;

    // Tests
    setup.check_cp = true;
    setup.check_iter = true;
    setup.check_stop = true;
    setup.check_safeguard_alpha = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
