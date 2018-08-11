__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

from Optizelle.Utility import *
from Optizelle.Properties import *
import Optizelle.InequalityConstrained.State

class X_Vectors(list):
    """Holds restart information for the vectors in the vector space X"""
    pass
class Z_Vectors(list):
    """Holds restart information for the vectors in the vector space Z"""
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

def release(X,Z,state,xs,zs,reals,nats,params):
    """Release the data into structures controlled by the user"""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Release the information from the state
    InequalityConstrainedRestartRelease(X,Z,state,xs,zs,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None

def capture(X,Z,state,xs,zs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    Optizelle.InequalityConstrained.State.checkT("state",state)
    checkVectors('xs',xs)
    checkVectors('zs',zs)
    checkReals('reals',reals)
    checkNaturals('nats',nats)
    checkParams('params',params)

    # Capture the restart information
    InequalityConstrainedRestartCapture(X,Z,state,xs,zs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None
