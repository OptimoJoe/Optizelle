// Test the quasinormal step where we cut off the Newton step with both a
// trust-region and a safeguard to force a dogleg step on the intersection of
// two circles constraint.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "qn.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess 
    auto x = std::vector <Real> { 0., 0. };
    auto y = std::vector <Real> { 0., 0. };

    // Setup the test 
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::CircleIntersection(1.,0.,1.,1.));
    
    // Unrestricted Cauchy point is (0.2,0.2) and Newton is (0,0.5)
    auto cp = std::vector <Real> {0.2,0.2};
    auto n = std::vector <Real> {0.0,0.5};

    // Set how far we want to traverse between them.  The first iteration sets
    // up the trust-region whereas the second sets up the inequality.
    auto dog = cp;
    for(unsigned i=0;i<=1;i++) {
        auto alpha = i==0 ? 0.7 : 0.6;
        X::copy(cp,dog);
        X::scal(Real(1.)-alpha,dog);
        X::axpy(alpha,n,dog);
        auto c = setup.zeta * setup.gamma;
        if(i==0)
            setup.delta = std::sqrt(X::innr(dog,dog)) / c;
        else
            setup.h.reset(new Unit <Real>::Constraint::Box(
                {-1.,-1.},{Real(2.),dog[1]/c}));
    }
    
    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::DoglegSafeguard;
    setup.dx_ncp_star.reset(new std::vector <Real> (cp));
    setup.dx_n_star.reset(new std::vector <Real> (dog)); 
    
    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;
    setup.check_safe = true;
    setup.check_augsys = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success 
    return EXIT_SUCCESS;
}
