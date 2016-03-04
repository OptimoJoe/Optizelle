// Forces the quasinormal step to exit on the first step where we're stuck
// on a local min for the Cauchy point formulation

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "qn.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess 
    auto x = std::vector <Real> { 0.0, 0.0 };
    auto y = std::vector <Real> { 0. };

    // Setup the test 
    auto setup = Unit <Real>::QN(x,y);
    setup.g.reset(new Unit <Real>::Constraint::Quadratic(1.,1.));

    // Set the targets
    setup.qn_stop_star = Optizelle::QuasinormalStop::LocalMin;
    setup.dx_n_star.reset(new std::vector <Real> {0.0,0.0});
    setup.dx_ncp_star = std::make_unique <X_Vector> (*setup.dx_n_star);
    
    // Set what tests we want
    setup.check_stop = true;
    setup.check_dx_n = true;
    setup.check_dx_ncp = true;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success 
    return EXIT_SUCCESS;
}
