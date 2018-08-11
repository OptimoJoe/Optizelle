__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization"

import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Functions import *

def getMin(X, Z, msg, fns, state, smanip=None):
    """Solves an inequality constrained optimization problem
    Basic solve: getMin(X,Z,msg,fns,state)
    Solve with a state manipulator: getMin(X,Z,msg,fns,state,smanip)
    """
    if smanip is None:
        smanip = StateManipulator()

    # Check the arguments
    checkVectorSpace("X",X)
    checkVectorSpace("Z",Z)
    checkMessaging("msg",msg)
    Optizelle.InequalityConstrained.Functions.checkT("fns",fns)
    Optizelle.InequalityConstrained.State.checkT("state",state)
    checkStateManipulator("smanip",smanip)

    # Call the optimization
    InequalityConstrainedAlgorithmsGetMin(X,Z,msg,fns,state,smanip)
