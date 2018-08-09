__all__ = [
    "read",
    "write_restart",
    "read_restart"
]

import Optizelle.EqualityConstrained.State
from Optizelle.Utility import *
from Optizelle.Properties import *

def read(X,Y,fname,state):
    """Read parameters from file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the read
    Optizelle.Utility.EqualityConstrainedStateReadJson(X,Y,fname,state)

def write_restart(X,Y,fname,state):
    """Writes a json restart file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the write
    EqualityConstrainedRestartWriteRestart(X,Y,fname,state)

def read_restart(X,Y,fname,x,y,state):
    """Reads a json restart file"""

    # Check our arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the read
    EqualityConstrainedRestartReadRestart(X,Y,fname,x,y,state)
