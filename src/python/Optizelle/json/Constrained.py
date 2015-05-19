__all__ = [
    "read",
    "write_restart",
    "read_restart"
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

def write_restart(X,Y,Z,msg,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the write
    Optizelle.Utility.ConstrainedRestartWriteRestart(
        X,Y,Z,msg,fname,state)

def read_restart(X,Y,Z,msg,fname,x,y,z,state):
    """Reads a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.Constrained.State.checkT("state",state)

    # Do the read 
    Optizelle.Utility.ConstrainedRestartReadRestart(
        X,Y,Z,msg,fname,x,y,z,state)
