__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.Constrained.State
import Optizelle.Constrained.Functions

def getMin(X, Y, Z, msg, fns, state, smanip=None):
    """Solves a constrained optimization problem
    Basic solve: getMin(X,Y,Z,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Y,Z,msg,smanip,fns,state)
    """
    if smanip is None:
        smanip = Optizelle.StateManipulator()

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Y",Y)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.Constrained.Functions.checkT("fns",fns)
    Optizelle.Constrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.ConstrainedAlgorithmsGetMin(X,Y,Z,msg,fns,state,smanip)
