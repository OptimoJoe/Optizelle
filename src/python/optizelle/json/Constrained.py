__all__ = [
    "read"
]

import Optizelle
import Optizelle.Utility
import Optizelle.Constrained.State

def read(X,Y,Z,msg,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the read
    Optizelle.Utility.ConstrainedStateReadJson(X,Y,Z,msg,fname,state)
