// Test the quasinormal step where we cut off the Newton step and force a
// dogleg step on the intersection of two circles constraint.  In addition, we
// cut off the Cauchy point and make sure that it's infeasible.

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
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

    // Set how far we want to traverse between them
    auto alpha = 0.6;
    auto dog = cp;
    X::scal(Real(1.)-alpha,dog);
    X::axpy(alpha,n,dog);
    setup.delta = std::sqrt(X::innr(dog,dog)) / setup.zeta;

    // Create an inequality constraint to cut off the Cauchy point 
    auto c = setup.zeta * setup.gamma;
    setup.h.reset(new Unit <Real>::Constraint::Box(
        {-1.,-1.},{0.1/c,10.}));

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::DoglegTrustRegion;
    setup.dx_ncp_star.reset(new std::vector <Real> ({0.1,0.1}));
    setup.dx_n_star.reset(new std::vector <Real> (dog)); 
    
    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;
    setup.check_tr = true;
    setup.check_augsys = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success 
    return EXIT_SUCCESS;
}
