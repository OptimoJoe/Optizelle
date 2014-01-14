__all__ = [
    "read",
    "write_restart",
    "read_restart"
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

def write_restart(X,Y,msg,fname,state):
    """Writes a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the write
    Optizelle.Utility.EqualityConstrainedRestartWriteRestart(
        X,Y,msg,fname,state)

def read_restart(X,Y,msg,fname,x,y,state):
    """Reads a json restart file"""
    
    # Check our arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.checkString("fname",fname)
    Optizelle.EqualityConstrained.State.checkT("state",state)

    # Do the read 
    Optizelle.Utility.EqualityConstrainedRestartReadRestart(
        X,Y,msg,fname,x,y,state)
