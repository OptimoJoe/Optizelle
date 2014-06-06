__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.Constrained.State
import Optizelle.Constrained.Functions

def getMin(*args):
    """Solves a constrained optimization problem
    Basic solve: getMin(X,Y,Z,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Y,Z,msg,smanip,fns,state)
    """

    # Check the number of arguments
    if len(args)!=6 and len(args)!=7:
        raise Exception("The getMin function requires either 6 or 7 arguments, "
            "but %d given." % len(args))

    # Extract the arguments
    X=args[0]
    Y=args[1]
    Z=args[2]
    msg=args[3]
    fns = args[4] 
    state = args[5] 
    smanip = Optizelle.StateManipulator() if len(args)==6 else args[6]

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
