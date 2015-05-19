__all__ = [
    "t"
]

import Optizelle.Unconstrained.Functions
from Optizelle import \
    checkDelete,\
    checkScalarValuedFunction, \
    checkVectorValuedFunction, \
    checkOperator

class t(Optizelle.Unconstrained.Functions.t):
    """All the functions required by an optimization algorithm""" 
    def __init__(self):
        super(t,self).__init__()
        self._g=Optizelle.VectorValuedFunction()
        self._PSchur_left=Optizelle.Operator()
        self._PSchur_right=Optizelle.Operator()
    
    # Create all of the properties
    g = Optizelle.createVectorValuedFunctionProperty(
        "g",
        "Equality constraints")
    PSchur_left = Optizelle.createOperatorProperty(
        "PSchur_left",
        "Left preconditioner for the augmented system")
    PSchur_right = Optizelle.createOperatorProperty(
        "PSchur_right",
        "Right preconditioner for the augmented system")

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type EqualityConstrained.Functions.t."
            % (name))
