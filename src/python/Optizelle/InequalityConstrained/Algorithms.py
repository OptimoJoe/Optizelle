__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions

def getMin(X, Z, msg, fns, state, smanip=None):
    """Solves an inequality constrained optimization problem
    Basic solve: getMin(X,Z,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Z,msg,fns,state,smanip)
    """
    if smanip is None:
        smanip = Optizelle.StateManipulator()

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.InequalityConstrained.Functions.checkT("fns",fns)
    Optizelle.InequalityConstrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.InequalityConstrainedAlgorithmsGetMin(
        X,Z,msg,fns,state,smanip)
