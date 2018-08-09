__all__ = [
    "t"
]

from Optizelle.Properties import *
from Optizelle.Functions import *

class t(object):
    """All the functions required by an optimization algorithm"""
    def __init__(self):
        self._f=ScalarValuedFunction()
        self._PH=Operator()

    # Create all of the properties
    f = createScalarValuedFunctionProperty(
        "f",
        "Objective function")
    PH = createOperatorProperty(
        "PH",
        "Preconditioner for the Hessian of the objective")

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Unconstrained.Functions.t."
            % (name))
