__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Functions import *

def getMin(X, Y, msg, fns, state, smanip=None):
    """Solves an equality constrained optimization problem
    Basic solve: getMin(X,Y,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Y,msg,fns,state,smanip)
    """

    if smanip is None:
        smanip = Optizelle.StateManipulator()

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Y",Y)
    checkMessaging("msg",msg)
    Optizelle.EqualityConstrained.Functions.checkT("fns",fns)
    Optizelle.EqualityConstrained.State.checkT("state",state)
    checkStateManipulator("smanip",smanip)

    # Call the optimization
    EqualityConstrainedAlgorithmsGetMin(X,Y,msg,fns,state,smanip)
