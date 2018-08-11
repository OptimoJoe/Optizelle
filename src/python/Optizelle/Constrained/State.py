__all__ = [
    "t"
]

import Optizelle.EqualityConstrained.State
import Optizelle.InequalityConstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Enumerated import *

class t(
    Optizelle.EqualityConstrained.State.t,
    Optizelle.InequalityConstrained.State.t):
    """Internal state of the optimization"""

    def __init__(self,X,Y,Z,x,y,z):
        """Constructor"""

        # Check our arguments
        checkVectorSpace("X",X)
        checkVectorSpace("Y",Y)
        checkEuclidean("Z",Z)

        # Allocate memory for our vectors
        Optizelle.Unconstrained.State.allocateVectors(self,X,x)
        Optizelle.EqualityConstrained.State.allocateVectors(self,X,Y,x,y)
        Optizelle.InequalityConstrained.State.allocateVectors(self,X,Z,x,z)

        # Create the state
        ConstrainedStateCreate(self,X,Y,Z,x,y,z)

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Constrained.State.t."
            % (name))
