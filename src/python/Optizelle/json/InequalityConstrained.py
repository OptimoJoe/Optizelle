__all__ = [
    "read",
    "write_restart",
    "read_restart"
]

import Optizelle.Unconstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *

def read(X,Z,fname,state):
    """Read parameters from file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the read
    InequalityConstrainedStateReadJson(X,Z,fname,state)

def write_restart(X,Z,fname,state):
    """Writes a json restart file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the write
    InequalityConstrainedRestartWriteRestart(X,Z,fname,state)

def read_restart(X,Z,fname,x,z,state):
    """Reads a json restart file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the read
    InequalityConstrainedRestartReadRestart(X,Z,fname,x,z,state)
