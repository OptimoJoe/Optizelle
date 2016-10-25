__all__ = [
    "read",
    "write_restart",
    "read_restart"
]

import Optizelle.Unconstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *

def read(X,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    checkVectorSpace("X",X)
    checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the read
    UnconstrainedStateReadJson(X,fname,state)

def write_restart(X,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    checkVectorSpace("X",X)
    checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the write
    UnconstrainedRestartWriteRestart(X,fname,state)

def read_restart(X,fname,x,state):
    """Reads a json restart file"""
    
    # Check our arguments
    checkVectorSpace("X",X)
    checkString("fname",fname)
    Optizelle.Unconstrained.State.checkT("state",state)

    # Do the read 
    UnconstrainedRestartReadRestart(X,fname,x,state)
