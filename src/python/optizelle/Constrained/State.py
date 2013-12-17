__all__ = [
    "t"
]

import Optizelle.EqualityConstrained.State
import Optizelle.InequalityConstrained.State

class t(
    Optizelle.EqualityConstrained.State.t,
    Optizelle.InequalityConstrained.State.t):
    """Internal state of the optimization"""

    def __init__(self,X,Y,Z,msg,x,y,z):
        """Constructor"""

        # Check our arguments
        Optizelle.checkVectorSpace("X",X)
        Optizelle.checkVectorSpace("Y",Y)
        Optizelle.checkEuclidean("Z",Z)
        Optizelle.checkMessaging("msg",msg)

        # Create the state
        Optizelle.Utility.ConstrainedStateCreate(self,X,Y,Z,msg,x,y,z)

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Constrained.State.t."
            % (name))
