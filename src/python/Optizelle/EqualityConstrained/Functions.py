__all__ = [
    "t"
]

import Optizelle.Unconstrained.Functions
from Optizelle.Properties import *
from Optizelle.Functions import *

class t(Optizelle.Unconstrained.Functions.t):
    """All the functions required by an optimization algorithm"""
    def __init__(self):
        super(t,self).__init__()
        self._g=VectorValuedFunction()
        self._PSchur_left=Operator()
        self._PSchur_right=Operator()

    # Create all of the properties
    g = createVectorValuedFunctionProperty(
        "g",
        "Equality constraints")
    PSchur_left = createOperatorProperty(
        "PSchur_left",
        "Left preconditioner for the augmented system")
    PSchur_right = createOperatorProperty(
        "PSchur_right",
        "Right preconditioner for the augmented system")

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type EqualityConstrained.Functions.t."
            % (name))
