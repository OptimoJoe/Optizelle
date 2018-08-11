#include "optizelle/vspaces.h"

namespace Optizelle {
    // Different cones used in SQL problems
    namespace Cone {

        // Converts the cone to a string
        std::string to_string(t const & cone){
            switch(cone){
            case Linear:
                return "Linear";
            case Quadratic:
                return "Quadratic";
            case Semidefinite:
                return "Semidefinite";
            default:
                throw;
            }
        }

        // Converts a string to a cone
        t from_string(std::string const & cone){
            if(cone=="Linear")
                return Linear;
            else if(cone=="Quadratic")
                return Quadratic;
            else if(cone=="Semidefinite")
                return Semidefinite;
            else
                throw;
        }

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name) {
            if( name=="Linear" ||
                name=="Quadratic" ||
                name=="Semidefinite"
            )
                return true;
            else
                return false;
        }
    }

    // Optimization problems instantiated on these vector spaces.  In theory,
    // this should help our compilation times.
    template struct Unconstrained<double,Rm>;
    template struct Unconstrained<float,Rm>;
    template struct EqualityConstrained<double,Rm,Rm>;
    template struct EqualityConstrained<float,Rm,Rm>;
    template struct InequalityConstrained<double,Rm,Rm>;
    template struct InequalityConstrained<float,Rm,Rm>;
    template struct InequalityConstrained<double,Rm,SQL>;
    template struct InequalityConstrained<float,Rm,SQL>;
    template struct Constrained<double,Rm,Rm,Rm>;
    template struct Constrained<float,Rm,Rm,Rm>;
    template struct Constrained<double,Rm,Rm,SQL>;
    template struct Constrained<float,Rm,Rm,SQL>;
}
