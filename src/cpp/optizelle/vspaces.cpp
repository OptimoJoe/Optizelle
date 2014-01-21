/*
Copyright 2013 OptimoJoe.

For the full copyright notice, see LICENSE.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Joseph Young (joe@optimojoe.com)
*/

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
