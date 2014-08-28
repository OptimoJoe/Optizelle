__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions

def getMin(X, Y, msg, fns, state, smanip=None):
    """Solves an equality constrained optimization problem
    Basic solve: getMin(X,Y,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Y,msg,fns,state,smanip)
    """

    if smanip is None:
        smanip = Optizelle.StateManipulator()

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.EqualityConstrained.Functions.checkT("fns",fns)
    Optizelle.EqualityConstrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.EqualityConstrainedAlgorithmsGetMin(
        X,Y,msg,fns,state,smanip)
