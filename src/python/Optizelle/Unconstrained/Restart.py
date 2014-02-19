__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

import Optizelle
import Optizelle.Utility
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
    Optizelle.checkVectorSpace("X",X)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Release the information from the state
    Optizelle.Utility.UnconstrainedRestartRelease(
        X,Optizelle.Messaging(),state,xs,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None 

def capture(X,msg,state,xs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.Unconstrained.State.checkT("state",state)
    Optizelle.checkVectors('xs',xs)
    Optizelle.checkReals('reals',reals)
    Optizelle.checkNaturals('nats',nats)
    Optizelle.checkParams('params',params)

    # Capture the restart information
    Optizelle.Utility.UnconstrainedRestartCapture(
        X,msg,state,xs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None 
