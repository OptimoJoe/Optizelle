__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

from Optizelle.Utility import *
from Optizelle.Properties import *
import Optizelle.Constrained.State

class X_Vectors(list):
    """Holds restart information for the vectors in the vector space X"""
    pass
class Y_Vectors(list):
    """Holds restart information for the vectors in the vector space Y"""
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

def release(X,Y,Z,state,xs,ys,zs,reals,nats,params):
    """Release the data into structures controlled by the user"""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    Optizelle.Constrained.State.checkT("state",state)

    # Release the information from the state
    ConstrainedRestartRelease(X,Y,Z,state,xs,ys,zs,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None

def capture(X,Y,Z,state,xs,ys,zs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    Optizelle.Constrained.State.checkT("state",state)
    checkVectors('xs',xs)
    checkVectors('ys',ys)
    checkVectors('zs',zs)
    checkReals('reals',reals)
    checkNaturals('nats',nats)
    checkParams('params',params)

    # Capture the restart information
    ConstrainedRestartCapture(X,Y,Z,state,xs,ys,zs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None
