#include "peopt/vspaces.h"

namespace peopt {
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
