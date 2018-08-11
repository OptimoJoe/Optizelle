# Tests our ability to throw a native exception and catch it

import Optizelle
import numpy

# Define a function that just throws an exception
class Objective(Optizelle.ScalarValuedFunction):
    # Evaluation
    def eval(self,x):
        raise RuntimeError("Evaluation")

    # Gradient
    def grad(self,x,grad):
        raise RuntimeError("Gradient")

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        raise RuntimeError("Hessian-vector product")

# Allocate memory for an initial guess
x = numpy.array([1.2,2.3])

# Create an optimization state
state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x)

# Create a bundle of functions
fns=Optizelle.Unconstrained.Functions.t()
fns.f=Objective()

# Try to catch the error
try:
    Optizelle.Unconstrained.Algorithms.getMin(
        Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)
except RuntimeError as e:
    if e.message!="Gradient":
        raise
