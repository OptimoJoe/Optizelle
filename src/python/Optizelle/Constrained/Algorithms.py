__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle.Constrained.State
import Optizelle.Constrained.Functions
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Functions import *

def getMin(X, Y, Z, msg, fns, state, smanip=None):
    """Solves a constrained optimization problem
    Basic solve: getMin(X,Y,Z,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Y,Z,msg,smanip,fns,state)
    """
    if smanip is None:
        smanip = StateManipulator()

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkVectorSpace("Z",Z)
    checkMessaging("msg",msg)
    Optizelle.Constrained.Functions.checkT("fns",fns)
    Optizelle.Constrained.State.checkT("state",state)
    checkStateManipulator("smanip",smanip)

    # Call the optimization
    ConstrainedAlgorithmsGetMin(X,Y,Z,msg,fns,state,smanip)
