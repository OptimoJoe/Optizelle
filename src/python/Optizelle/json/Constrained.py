__all__ = [
    "read",
    "write_restart",
    "read_restart"
]

import Optizelle.Unconstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *

def read(X,Y,Z,fname,state):
    """Read parameters from file"""
        
    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the read
    ConstrainedStateReadJson(X,Y,Z,fname,state)

def write_restart(X,Y,Z,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the write
    ConstrainedRestartWriteRestart(X,Y,Z,fname,state)

def read_restart(X,Y,Z,fname,x,y,z,state):
    """Reads a json restart file"""
    
    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the read 
    ConstrainedRestartReadRestart(X,Y,Z,fname,x,y,z,state)
