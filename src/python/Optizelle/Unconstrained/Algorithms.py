__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization"

import Optizelle.Unconstrained.State
import Optizelle.Unconstrained.Functions
from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Functions import *

def getMin(X, msg, fns, state, smanip=None):
    """Solves an unconstrained optimization problem
    Basic solve: getMin(X,msg,fns,state)
    Solve with a state manipulator: getMin(X,msg,fns,state,smanip)
    """
    if smanip is None:
        smanip = StateManipulator()

    # Check the arguments
    checkVectorSpace("X",X)
    checkMessaging("msg",msg)
    Optizelle.Unconstrained.Functions.checkT("fns",fns)
    Optizelle.Unconstrained.State.checkT("state",state)
    checkStateManipulator("smanip",smanip)

    # Call the optimization
    UnconstrainedAlgorithmsGetMin(X,msg,fns,state,smanip)
