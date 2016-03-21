// Tests the nullspace projection when provided with a zero direction 

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "augsys.h"
#include "spaces.h"

int main(int argc,char* argv[]){

    // Generate an initial guess 
    auto x = std::vector <Real> { 0., 0.};
    auto y = std::vector <Real> { 0. };

    // Setup the test 
    auto setup = Unit <Real>::Proj(x,y);
    setup.g.reset(new Unit <Real>::Constraint::Linear);
    setup.dx.reset(new std::vector <Real> {1e-16,1e-16});

    // Set the targets
    setup.P_dx_star.reset(new std::vector <Real> {0.0,0.0});
    
    // Set what tests we want
    setup.check_sol = true;
    setup.check_null = true;
    setup.check_augsys = false;

    // Run the test
    Unit <Real>::run_and_verify(setup);

    // Declare success 
    return EXIT_SUCCESS;
}
