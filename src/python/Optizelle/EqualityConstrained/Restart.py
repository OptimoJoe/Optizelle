__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

from Optizelle.Utility import *
from Optizelle.Properties import *
import Optizelle.EqualityConstrained.State

class X_Vectors(list):
    """Holds restart information for the vectors in the vector space X"""
    pass
class Y_Vectors(list):
    """Holds restart information for the vectors in the vector space Y"""
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

def release(X,Y,state,xs,ys,reals,nats,params):
    """Release the data into structures controlled by the user"""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Release the information from the state
    EqualityConstrainedRestartRelease(X,Y,state,xs,ys,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None 

def capture(X,Y,state,xs,ys,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    Optizelle.EqualityConstrained.State.checkT("state",state)
    checkVectors('xs',xs)
    checkVectors('ys',ys)
    checkReals('reals',reals)
    checkNaturals('nats',nats)
    checkParams('params',params)

    # Capture the restart information
    EqualityConstrainedRestartCapture(X,Y,state,xs,ys,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None 
