__all__ = [
    "read"
]

import Optizelle
import Optizelle.Utility
import Optizelle.Unconstrained.State

def read(X,msg,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the read
    Optizelle.Utility.UnconstrainedStateReadJson(X,msg,fname,state)
