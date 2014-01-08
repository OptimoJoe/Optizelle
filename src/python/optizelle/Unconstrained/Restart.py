__all__ = [
    "release"
]

__doc__ = "Utilities for restarting the optimization"

import Optizelle
import Optizelle.Utility
import Optizelle.Unconstrained.State

def release(X,state):
    """Release the data into structures controlled by the user"""

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Create structures to hold the released information
    xs=([],[])
    reals=([],[])
    nats=([],[])
    params=([],[])

    # Release the information from the state
    Optizelle.Utility.UnconstrainedRestartRelease(
        X,Optizelle.Messaging(),state,xs,reals,nats,params)

    # Return the released information
    return (xs,reals,nats,params)

def capture(X,msg,state,xs,reals,nats,params):
    """Capture data from structures controlled by the user."""

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.Unconstrained.State.checkT("state",state)

    map(lambda x:Optizelle.checkString("xs[0] member",x),xs[0])
    map(lambda x:Optizelle.checkString("reals[0] member",x),reals[0])
    map(lambda x:Optizelle.checkString("nats[0] member",x),nats[0])
    map(lambda x:Optizelle.checkString("params[0] member",x),params[0])

    map(lambda x:Optizelle.checkFloat("reals[1] member",x),reals[1])
    map(lambda x:Optizelle.checkNatural("nats[1] member",x),nats[1])
    map(lambda x:Optizelle.checkString("params[1] member",x),params[1])

    # Capture the restart information
    Optizelle.Utility.UnconstrainedRestartCapture(
        X,msg,state,xs,reals,nats,params)

    # Return nothing.  The state has been modified.
    return None 
