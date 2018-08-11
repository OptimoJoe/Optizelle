__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

from Optizelle.Utility import *
from Optizelle.Properties import *
import Optizelle.Unconstrained.State

class X_Vectors(list):
    """Holds restart information for the vectors in the vector space X"""
    pass
class Reals(list):
    """Holds restart information for real numbers"""
    pass
class Naturals(list):
    """Holds restart information for natural numbers"""
    pass
class Params(list):
    """Holds restart information for parameters"""
    pass

def release(X,state,xs,reals,nats,params):
    """Release the data into structures controlled by the user"""

    # Check the arguments
    checkVectorSpace("X",X)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Release the information from the state
    UnconstrainedRestartRelease(X,state,xs,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None

def capture(X,state,xs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    checkVectorSpace("X",X)
    Optizelle.Unconstrained.State.checkT("state",state)
    checkVectors('xs',xs)
    checkReals('reals',reals)
    checkNaturals('nats',nats)
    checkParams('params',params)

    # Capture the restart information
    UnconstrainedRestartCapture(X,state,xs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None
