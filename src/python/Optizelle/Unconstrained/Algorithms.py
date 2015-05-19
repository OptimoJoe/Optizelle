__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.Unconstrained.State
import Optizelle.Unconstrained.Functions

def getMin(X, msg, fns, state, smanip=None):
    """Solves an unconstrained optimization problem
    Basic solve: getMin(X,msg,fns,state) 
    Solve with a state manipulator: getMin(X,msg,fns,state,smanip)
    """
    if smanip is None:
        smanip = Optizelle.StateManipulator()

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.Unconstrained.Functions.checkT("fns",fns)
    Optizelle.Unconstrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.UnconstrainedAlgorithmsGetMin(X,msg,fns,state,smanip)
