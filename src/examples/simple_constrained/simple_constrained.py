# Optimize a simple optimization problem with an optimal solution
# of (1/3,1/3)

import Optizelle 
import numpy
import sys

# Squares its input
sq = lambda x:x*x

# Define a simple objective where 
# 
# f(x,y)=(x+1)^2+(y+1)^2
#
class MyObj(Optizelle.ScalarValuedFunction):

    # Evaluation 
    def eval(self,x):
        return sq(x[0]+1.)+sq(x[1]+1.)

    # Gradient
    def grad(self,x,grad):
        grad[0]=2.*x[0]+2.
        grad[1]=2.*x[1]+2.

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        H_dx[0]=2.*dx[0]
        H_dx[1]=2.*dx[1]

# Define a simple equality
#
# g(x,y)= [ x + 2y = 1 ] 
#
class MyEq(Optizelle.VectorValuedFunction):

    # y=g(x) 
    def eval(self,x,y):
        y[0]=x[0]+2.*x[1]-1.

    # y=g'(x)dx
    def p(self,x,dx,y):
        y[0]= dx[0]+2.*dx[1]

    # xhat=g'(x)*dy
    def ps(self,x,dy,xhat):
        xhat[0]= dy[0]
        xhat[1]= 2.*dy[0]

    # xhat=(g''(x)dx)*dy
    def pps(self,x,dx,dy,xhat):
        xhat.fill(0.)

# Define simple inequalities 
#
# h(x,y)= [ 2x + y >= 1 ] 
#
class MyIneq(Optizelle.VectorValuedFunction):

    # z=h(x) 
    def eval(self,x,z):
        z[0]=2.*x[0]+x[1]-1.

    # z=h'(x)dx
    def p(self,x,dx,z):
        z[0]= 2.*dx[0]+dx[1]

    # xhat=h'(x)*dz
    def ps(self,x,dz,xhat):
        xhat[0]= 2.*dz[0]
        xhat[1]= dz[0]

    # xhat=(h''(x)dx)*dz
    def pps(self,x,dx,dz,xhat):
        xhat.fill(0.)

# Read in the name for the input file
if len(sys.argv) != 2:
    sys.stderr.write("Usage: python simple_constrained.py <parameters>.json\n")
    raise ValueError("Parameters JSON file required.")

fname = sys.argv[1]

# Generate an initial guess 
x = numpy.array([2.1,1.1])

# Allocate memory for the equality multiplier 
y = numpy.array([0.])

# Allocate memory for the inequality multiplier 
z = numpy.array([0.])

# Create an optimization state
state=Optizelle.Constrained.State.t(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z)

# Read the parameters from file
Optizelle.json.Constrained.read(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,fname,state)

# Create a bundle of functions
fns=Optizelle.Constrained.Functions.t()
fns.f=MyObj()
fns.g=MyEq()
fns.h=MyIneq()

# Solve the optimization problem
Optizelle.Constrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])

# Write out the final answer to file
Optizelle.json.Constrained.write_restart(Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,
    "solution.json",state)
