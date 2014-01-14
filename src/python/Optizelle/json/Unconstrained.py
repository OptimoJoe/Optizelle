__all__ = [
    "read",
    "write_restart",
    "read_restart"
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

def write_restart(X,msg,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the write
    Optizelle.Utility.UnconstrainedRestartWriteRestart(X,msg,fname,state)

def read_restart(X,msg,fname,x,state):
    """Reads a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the read 
    Optizelle.Utility.UnconstrainedRestartReadRestart(X,msg,fname,x,state)
