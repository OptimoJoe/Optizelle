__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions

def getMin(*args):
    """Solves an inequality constrained optimization problem
    Basic solve: getMin(X,Z,msg,fns,state) 
    Solve with a state manipulator: getMin(X,Z,msg,smanip,fns,state)
    """

    # Check the number of arguments
    if len(args)!=5 and len(args)!=6:
        raise Exception("The getMin function requires either 5 or 6 arguments, "
            "but %d given." % len(args))

    # Extract the arguments
    X=args[0]
    Z=args[1]
    msg=args[2]
    fns = args[3] if len(args)==5 else args[4]
    state = args[4] if len(args)==5 else args[5]
    smanip = Optizelle.StateManipulator() if len(args)==5 else args[3]

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkVectorSpace("Z",Z)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.InequalityConstrained.Functions.checkT("fns",fns)
    Optizelle.InequalityConstrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.InequalityConstrainedAlgorithmsGetMin(
        X,Z,msg,smanip,fns,state)
