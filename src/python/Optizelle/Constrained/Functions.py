__all__ = [
    "t"
]

import Optizelle.EqualityConstrained.Functions
import Optizelle.InequalityConstrained.Functions

class t(
    Optizelle.EqualityConstrained.Functions.t,
    Optizelle.InequalityConstrained.Functions.t):
    """All the functions required by an optimization algorithm"""
    def __init__(self):
        super(t,self).__init__()

def checkT(name,value):
    """Check that we have a bundle of functions"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Constrained.Functions.t."
            % (name))
