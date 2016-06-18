//  Tests error handling when reading invalid parameters 

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "optizelle/json.h"
    
// Create some type shortcuts
template <typename Real> using XX = Optizelle::Rm <Real>;
typedef double Real;

int main() {
    // Set the parameter file name 
    auto fname = "bad_params.json"; 

    // Create a type shortcut
    using Optizelle::Rm;

    // Allocate memory for an initial guess
    auto x = std::vector <double> {1.2, 2.3};
    
    // Create an optimization state
    Optizelle::Unconstrained<double,Rm>::State::t state(x);

    auto msg = std::string("");
    //---Exception0---
    // Read parameters from file
    try {
        Optizelle::json::Unconstrained <Real,XX>::read(fname,state);
    } catch(Optizelle::Exception::t const & e) {
        // Convert the error message to a string 
        msg = Optizelle::Exception::to_string(e);

        // Print the error message directly 
        Optizelle::Exception::to_stderr(e);
    }
    //---Exception1---
    
    // If we didn't throw an exception above, throw an error 
    if(msg.length()==0)
        throw Optizelle::Exception::t(
            "Error catching missed our bad parameter");

    return EXIT_SUCCESS;
}
