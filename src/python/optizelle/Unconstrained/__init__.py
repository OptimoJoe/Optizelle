import Optizelle

__all__ = [
    "State",
    "Functions",
    "Algorithms"
]
__doc__ = "Unconstrained optimization problems."

import Optizelle

class X_Vectors(Optizelle.RestartPackage):
    """Holds restart information for the vectors in the vector space X"""
    pass
class reals(Optizelle.RestartPackage):
    """Holds restart information for real numbers"""
    pass
class nats(Optizelle.RestartPackage):
    """Holds restart information for natural numbers"""
    pass
class params(Optizelle.RestartPackage):
    """Holds restart information for parameters"""
    pass
