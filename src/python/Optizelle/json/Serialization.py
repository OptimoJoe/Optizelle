__all__ = [
    "serialize",
    "deserialize"
]
__doc__ = "Interaction between Optizelle with JSON formatted files"

import numpy
import Optizelle.Exception
from Optizelle.Properties import *

class Extendable(object):
    """Allows a function to be extended"""
    def __init__(self, fn):
        self.fns = {}
        self.fn = fn

    def register(self,fn,vector_type):
        """Extends the current function with fn.  This function will only be 
        called if the first argument matches vector_type."""

        # Check our arguments
        checkFunction("(de)serialize",fn)
        checkType("vector_type",vector_type)

        # Register the function
        self.fns[vector_type]=fn

    def __call__(self,*args):
        try:
            return self.fns[type(args[0])](*args)
        except:
            return self.fn(*args) 

@ Extendable
def serialize(x,name,iter):
    """Converts a vector to a JSON formatted string"""
    
    raise Optizelle.Exception.t(
        "The serialize function for the vector %s not defined." % str(x))
    
@ Extendable
def deserialize(x,x_json):
    """Converts a JSON formatted string to a vector"""

    raise Optizelle.Exception.t(
        "The deserialize function for the vector %s not defined." % str(x))

def serialize_Rm(x,name,iter):
    """Serializes a numpy array for the vector space Optizelle.Rm""" 

    # Create the json representation
    x_json="[ "
    for i in xrange(x.size):
        x_json  += str(x[i]) + ", "
    x_json=x_json[0:-2]
    x_json +=" ]"

    return x_json

def deserialize_Rm(x,x_json):
    """Deserializes a numpy array for the vector space Optizelle.Rm""" 

    # Eliminate all whitespace
    x_json="".join(x_json.split())

    # Check if we're a vector
    if x_json[0:1]!="[" or x_json[-1:]!="]":
        raise TypeError("Attempted to deserialize a non-numpy.array vector.")

    # Eliminate the initial and final delimiters
    x_json=x_json[1:-1]

    # Create a list of the numbers involved 
    x_json=x_json.split(",")

    # Convert the strings to numbers
    x_json=map(lambda x:float(x),x_json)

    # Create an Optizelle.Rm vector
    return numpy.array(x_json)

# Register the serialization routines for numpy arrays 
serialize.register(serialize_Rm,numpy.ndarray)
deserialize.register(deserialize_Rm,numpy.ndarray)
