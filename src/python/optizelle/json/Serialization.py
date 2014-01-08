__all__ = [
    "serialize",
    "deserialize"
]
__doc__ = "Interaction between Optizelle with JSON formatted files"

import Optizelle
import numpy

class Extendable(object):
    """Allows a function to be extended"""
    def __init__(self, fn):
        self.fns = []
        self.fn = fn
    def register(self,fn):
        self.fns.append(fn)
    def __call__(self,x):
        for f in self.fns:
            try:
                return f(x)
            except:
                pass
        return self.fn(x)

@ Extendable
def serialize(x):
    """Converts a vector to a JSON formatted string"""
    raise Optizelle.Exception(
        "The serialize function for the vector %s not defined." % str(x))
    
@ Extendable
def deserialize(x):
    """Converts a JSON formatted string to a vector"""
    raise Optizelle.Exception(
        "The deserialize function for the vector %s not defined." % str(x))

def serialize_Rm(x):
    """Serializes a numpy array for the vector space Optizelle.Rm""" 

    # Check if we have a numpy array 
    if type(x)!=numpy.ndarray:
        raise TypeError("Attempted to serialize a non-numpy.array vector.")

    # Create the json representation
    x_json="{ [ "
    for i in xrange(x.size):
        x_json  += str(x[i]) + ", "
    x_json=x_json[0:-2]
    x_json +=" ] }"
    return x_json

def deserialize_Rm(x_json):
    """Deserializes a numpy array for the vector space Optizelle.Rm""" 

    # Eliminate all whitespace
    x_json="".join(x_json.split())

    # Check if we're a vector
    if x_json[0:2]!="{[" or x_json[-2:]!="]}":
        raise TypeError("Attempted to deserialize a non-numpy.array vector.")

    # Eliminate the initial and final delimiters
    x_json=x_json[2:-2]

    # Create a list of the numbers involved 
    x_json=x_json.split(",")

    # Convert the strings to numbers
    x_json=map(lambda x:float(x),x_json)

    # Create an Optizelle.Rm vector
    return numpy.array(x_json)

# Register the serialization routines for numpy arrays 
serialize.register(serialize_Rm)
deserialize.register(deserialize_Rm)
