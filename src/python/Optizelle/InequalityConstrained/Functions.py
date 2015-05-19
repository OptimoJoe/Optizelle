__all__ = [
    "t"
]

import Optizelle
import Optizelle.Unconstrained.Functions

class t(Optizelle.Unconstrained.Functions.t):
    """All the functions required by an optimization algorithm""" 
    def __init__(self):
        super(t,self).__init__()
        self._h=Optizelle.VectorValuedFunction()
    
    # Create all of the properties
    h = Optizelle.createVectorValuedFunctionProperty(
        "h",
        "Inequality constraints")

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type InequalityConstrained.Functions.t."
            % (name))
