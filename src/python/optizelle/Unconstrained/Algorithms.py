__all__ = [
    "getMin"
]

__doc__ = "Different algorithms used for optimization" 

import Optizelle
import Optizelle.Utility
import Optizelle.Unconstrained.State
import Optizelle.Unconstrained.Functions

def getMin(*args):
    """Solves an unconstrained optimization problem
    Basic solve: getMin(X,msg,fns,state) 
    Solve with a state manipulator: getMin(X,msg,smanip,fns,state)
    """

    # Check the number of arguments
    if len(args)!=4 and len(args)!=5:
        raise Exception("The getMin function requires either 4 or 5 arguments, "
            "but %d given." % len(args))

    # Extract the arguments
    X=args[0]
    msg=args[1]
    fns = args[2] if len(args)==4 else args[3]
    state = args[3] if len(args)==4 else args[4]
    smanip = Optizelle.StateManipulator() if len(args)==4 else args[2]

    # Check the arguments
    Optizelle.checkVectorSpace("X",X)
    Optizelle.checkMessaging("msg",msg)
    Optizelle.Unconstrained.Functions.checkT("fns",fns)
    Optizelle.Unconstrained.State.checkT("state",state)
    Optizelle.checkStateManipulator("smanip",smanip)

    # Call the optimization
    Optizelle.Utility.UnconstrainedAlgorithmsGetMin(X,msg,smanip,fns,state)
