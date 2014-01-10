__all__ = [
    "read",
    "write_restart",
    "read_restart"
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

def write_restart(X,Z,msg,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the write
    Optizelle.Utility.InequalityConstrainedRestartWriteRestart(
        X,Z,msg,fname,state)

def read_restart(X,Z,msg,fname,x,z,state):
    """Reads a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.InequalityConstrained.State.checkT("state",state)

    # Do the read 
    Optizelle.Utility.InequalityConstrainedRestartReadRestart(
        X,Z,msg,fname,x,z,state)
