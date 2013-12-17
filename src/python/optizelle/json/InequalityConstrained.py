__all__ = [
    "read"
]

import Optizelle
import Optizelle.Utility
import Optizelle.InequalityConstrained.State

def read(X,Z,msg,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the read
    Optizelle.Utility.InequalityConstrainedStateReadJson(X,Z,msg,fname,state)
