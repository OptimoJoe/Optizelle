// Run TCG with a safeguard where we have too many failed safeguard steps
// and are forced to both exit and retreat to a prior point.  Specifically
//
// x_offset = (0,0,0,0,0)
// w = (0,0,0,0,-1)
// lb = -7/12
//
// The sequence of solutions produced by CG is approximately
//
// 0 0 0 0 0
// 0 0 0 0 1/2
// 0 0 0 -1/3 2/3
// 0 0 1/4 -1/2 3/4
// 0 -1/5 2/5 -3/5 4/5
// 67000 -53000 40000 -27000 13000
//
// Hence, the sequence <x + x_offset,w> is
//
// -1/2, -2/3, -3/4, -4/5, -13000 
//
// As such, we require that <x + x_offset,w> >= -7/12 which means that we
// violate the safeguard afer the first iteration. 

#include "linear_algebra.h"
#include "spaces.h"

int main() {
    // Setup the problem 
    auto setup = Unit::tcg <Real,Rm> ();

    // Problem setup 
    setup.A = std::make_unique <Matrix>(
        Unit::Matrix <Real>::mostly_dd_indef(setup.m));
    setup.b = std::make_unique <Vector> (Unit::Vector <Real>::elast(setup.m));
    auto lb = Real(-7./12.); 
    auto x = Unit::Vector<Real>::zero(setup.m);
    auto w = std::vector <Real> {0,0,0,0,-1};
    setup.safeguard = std::make_unique<Optizelle::SafeguardSimplified<Real,Rm>>(
        Unit::Safeguard <Real,Rm>::lower(x,lb,w));
    setup.failed_max = 2;
    setup.delta = 1e5; 

    // Target solutions
    setup.iter_star = 3;
    setup.stop_star = Optizelle::TruncatedStop::TooManyFailedSafeguard;
    setup.safeguard_failed_star = 2;

    // Tests
    Unit::reset_checks <Real,Rm> (setup);
    setup.check_iter = true;
    setup.check_stop = true;
    setup.check_safeguard_failed = true;
    setup.check_safeguard_alpha = true;

    // Check the solver 
    Unit::run_and_verify <Real,Rm> (setup);

    // Declare success
    return EXIT_SUCCESS;
}
