__all__ = [
    "t"
]

import Optizelle.Unconstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Enumerated import *

def allocateVectors(self,X,Z,x,z):
    """Allocates memory for the state vectors"""
    self.z=Z.init(z)
    self.dz=Z.init(z)
    self.h_x=Z.init(z)

class t(Optizelle.Unconstrained.State.t):
    """Internal state of the optimization"""

    def __init__(self,X,Z,x,z):
        """Constructor"""

        # Check our arguments
        checkVectorSpace("X",X)
        checkEuclidean("Z",Z)
        
        # Allocate memory for our vectors
        Optizelle.Unconstrained.State.allocateVectors(self,X,x)
        allocateVectors(self,X,Z,x,z)

        # Create the state
        InequalityConstrainedStateCreate(self,X,Z,x,z)

    # Create all of the properties
    z = createVectorProperty(
        "z",
        "Inequality multiplier (dual variable or Lagrange multiplier)")
    dz = createVectorProperty(
        "dz",
        "Step in the inequality multiplier")
    h_x = createVectorProperty(
        "h_x",
        "The inequality constraint evaluated at x.")
    mu = createFloatProperty(
        "mu",
        "Interior point parameter")
    mu_est = createFloatProperty(
        "mu_est",
        "Current interior point estimate")
    mu_typ = createFloatProperty(
        "mu_typ",
        "Typical value for mu.  Generally, the first estimated value for mu.")
    eps_mu = createFloatProperty(
        "eps_mu",
        "Relative stopping criteria for the interior point parameter")
    sigma = createFloatProperty(
        "sigma",
        ("The amount that we reduce the interior point parameter by everytime "
        "we approach the central path"))
    gamma = createFloatProperty(
        "gamma",
        "How close we move to the boundary during a single step")
    alpha_z = createFloatProperty(
        "alpha_z",
        ("Amount we truncate dx in order to maintain feasibility "
        " of the inequality multiplier"))
    h_diag = createEnumProperty(
        "h_diag",
        FunctionDiagnostics,
        "Function diagnostics on h")
    z_diag = createEnumProperty(
        "z_diag",
        VectorSpaceDiagnostics,
        "Vector space diagnostics on Z")

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type InequalityConstrained.State.t."
            % (name))
