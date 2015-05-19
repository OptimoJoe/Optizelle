__all__ = [
    "t"
]

import Optizelle

class t(object):
    """All the functions required by an optimization algorithm""" 
    def __init__(self):
        self._f=Optizelle.ScalarValuedFunction()
        self._PH=Optizelle.Operator()

    # Create all of the properties
    f = Optizelle.createScalarValuedFunctionProperty(
        "f",
        "Objective function")
    PH = Optizelle.createOperatorProperty(
        "PH",
        "Preconditioner for the Hessian of the objective")

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Unconstrained.Functions.t."
            % (name))
