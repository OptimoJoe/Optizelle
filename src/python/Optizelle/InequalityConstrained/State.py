__all__ = [
    "t"
]

import Optizelle
import Optizelle.Unconstrained.State

def allocateVectors(self,X,Z,x,z):
    """Allocates memory for the state vectors"""
    self.z=Z.init(z)
    self.dz=Z.init(z)
    self.h_x=Z.init(z)

class t(Optizelle.Unconstrained.State.t):
    """Internal state of the optimization"""

    def __init__(self,X,Z,msg,x,z):
        """Constructor"""

        # Check our arguments
        Optizelle.checkVectorSpace("X",X)
        Optizelle.checkEuclidean("Z",Z)
        Optizelle.checkMessaging("msg",msg)
        
        # Allocate memory for our vectors
        Optizelle.Unconstrained.State.allocateVectors(self,X,x)
        allocateVectors(self,X,Z,x,z)

        # Create the state
        Optizelle.Utility.InequalityConstrainedStateCreate(self,X,Z,msg,x,z)

    # Create all of the properties
    z = Optizelle.createVectorProperty(
        "z",
        "Inequality multiplier (dual variable or Lagrange multiplier)")
    dz = Optizelle.createVectorProperty(
        "dz",
        "Step in the inequality multiplier")
    h_x = Optizelle.createVectorProperty(
        "h_x",
        "The inequality constraint evaluated at x.")
    mu = Optizelle.createFloatProperty(
        "mu",
        "Interior point parameter")
    mu_est = Optizelle.createFloatProperty(
        "mu_est",
        "Current interior point estimate")
    mu_typ = Optizelle.createFloatProperty(
        "mu_typ",
        "Typical value for mu.  Generally, the first estimated value for mu.")
    eps_mu = Optizelle.createFloatProperty(
        "eps_mu",
        "Relative stopping criteria for the interior point parameter")
    sigma = Optizelle.createFloatProperty(
        "sigma",
        ("The amount that we reduce the interior point parameter by everytime "
        "we approach the central path"))
    gamma = Optizelle.createFloatProperty(
        "gamma",
        "How close we move to the boundary during a single step")
    alpha_x = Optizelle.createFloatProperty(
        "alpha_x",
        ("Amount we truncate dx in order to maintain feasibility "
        " with respect to the inequality constraint"))
    alpha_z = Optizelle.createFloatProperty(
        "alpha_z",
        ("Amount we truncate dx in order to maintain feasibility "
        " of the inequality multiplier"))
    ipm = Optizelle.createEnumProperty(
        "ipm",
        Optizelle.InteriorPointMethod,
        "Type of interior point method")
    cstrat = Optizelle.createEnumProperty(
        "cstrat",
        Optizelle.CentralityStrategy,
        "Centrality strategy")
    h_diag = Optizelle.createEnumProperty(
        "h_diag",
        Optizelle.FunctionDiagnostics,
        "Function diagnostics on h")
    y_diag = Optizelle.createEnumProperty(
        "y_diag",
        Optizelle.VectorSpaceDiagnostics,
        "Vector space diagnostics on Y")
    delta_z = Optizelle.createFloatProperty(
        "delta_z",
        "Trust-region radius for the inequality constraint")

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type InequalityConstrained.State.t."
            % (name))
