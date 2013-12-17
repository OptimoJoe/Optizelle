__all__ = [
    "read"
]

import Optizelle
import Optizelle.Utility
import Optizelle.EqualityConstrained.State

def read(X,Y,msg,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the read
    Optizelle.Utility.EqualityConstrainedStateReadJson(X,Y,msg,fname,state)
