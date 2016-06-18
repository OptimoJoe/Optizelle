# A problem that helps us determine how to scale inequality constrained
# optimization problems

import Optizelle 
import numpy
import sys

# Define a simple objective where 
# 
# f(x) = 0.5 || x - c ||^2
#
class MyObj(Optizelle.ScalarValuedFunction):
    # Grab the center of the objective
    def __init__(self, c):
        self.c = c 

    # Evaluation 
    def eval(self,x):
        diff = x - self.c
        return 0.5*numpy.dot(diff,diff)

    # Gradient
    def grad(self,x,grad):
        numpy.copyto(grad,x-self.c);

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        numpy.copyto(H_dx,dx);

# Define simple inequalities 
#
# h(x) = x - lb 
#
class MyIneq(Optizelle.VectorValuedFunction):
    # Grab the lower bound 
    def __init__(self, lb):
        self.lb = lb 

    # z=h(x) 
    def eval(self,x,z):
        numpy.copyto(z,x-lb);

    # z=h'(x)dx
    def p(self,x,dx,z):
        numpy.copyto(z,dx);

    # xhat=h'(x)*dz
    def ps(self,x,dz,xhat):
        numpy.copyto(xhat,dz);

    # xhat=(h''(x)dx)*dz
    def pps(self,x,dx,dz,xhat):
        xhat.fill(0.)

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("inequality_scaling.py <parameters>")
fname=sys.argv[1]

# Set the size
m = 10;

# Generate an initial guess 
x = numpy.array([1.+10**(-x) for x in xrange(1,m+1)])

# Allocate memory for the inequality multiplier 
z = numpy.array(m*[0.])

# Create the center of the objective function
c = numpy.array(m*[-1.])

# Create the lower bound for the problem
lb = numpy.array(m*[1.])

# Create an optimization state
state=Optizelle.InequalityConstrained.State.t(Optizelle.Rm,Optizelle.Rm,x,z)

# Read the parameters from file
Optizelle.json.InequalityConstrained.read(
    Optizelle.Rm,Optizelle.Rm,fname,state)

# Create a bundle of functions
fns=Optizelle.InequalityConstrained.Functions.t()
fns.f=MyObj(c)
fns.h=MyIneq(lb)

# Solve the optimization problem
Optizelle.InequalityConstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: ["
for i in xrange(m):
    print "%1.16e" % state.x[i]
print "]"

# Write out the final answer to file
Optizelle.json.InequalityConstrained.write_restart(
    Optizelle.Rm,Optizelle.Rm,"solution.json",state)
