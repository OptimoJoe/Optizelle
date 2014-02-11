__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

import Optizelle
import Optizelle.Utility
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
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Release the information from the state
    Optizelle.Utility.InequalityConstrainedRestartRelease(
        X,Z,Optizelle.Messaging(),state,xs,zs,reals,nats,params)

    # Return nothing.  We've modified the passed in lists.
    return None 

def capture(X,Z,msg,state,xs,zs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.InequalityConstrained.State.checkT("state",state)
    Optizelle.checkVectors('xs',xs)
    Optizelle.checkVectors('zs',zs)
    Optizelle.checkReals('reals',reals)
    Optizelle.checkNaturals('nats',nats)
    Optizelle.checkParams('params',params)

    # Capture the restart information
    Optizelle.Utility.InequalityConstrainedRestartCapture(
        X,Z,msg,state,xs,zs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None 
