# Optimize a simple optimization problem with an optimal solution of (1/3,1/3),
#
# min x + y
# st  x + 2y >= 1
#     2x + y >= 1
#
# Now, in the case we don't have a starting feasible solution, we can play
# a reformulation trick that adds two scalar variables and allows us to find
# a strictly feasible solution.  Namely,
#
# min x + y
# st  x + 2y >= 1 - z
#     2x + y >= 1 - z
#     epsilon >= w             
#     z = w
#
# Note, most of the time, we're much better off just adding slack variables.
# Basically, this trick is only worthwhile when we don't have a linear system
# solver for the equality constraints added from the slacks since this method
# only adds a single equality constraint.

import Optizelle 
import numpy
import sys

# Squares its input
sq = lambda x:x*x

# Define an objective where 
# 
# f(x,y,z,w)=x+y
#
class MyObj(Optizelle.ScalarValuedFunction):
    # Evaluation 
    def eval(self,x):
        return x[0]+x[1];

    # Gradient
    def grad(self,x,grad):
        grad[0]=1.
        grad[1]=1.
        grad[2]=0.
        grad[3]=0.

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        H_dx.fill(0.)

# Define a single equality where
#
# g(x,y,z,w) = z - w = 0
#
class MyEq(Optizelle.VectorValuedFunction):

    # y=g(x) 
    def eval(self,x,y):
        y[0] = x[2]-x[3]

    # y=g'(x)dx
    def p(self,x,dx,y):
        y[0]= dx[2]-dx[3]

    # xhat=g'(x)*dy
    def ps(self,x,dy,xhat):
        xhat[0]= 0.
        xhat[1]= 0.
        xhat[2]= dy[0]
        xhat[3]= -dy[0]

    # xhat=(g''(x)dx)*dy
    def pps(self,x,dx,dy,xhat):
        xhat.fill(0.)

# Define some inequalities where
#
# h(x,y,z,w) = [ x + 2y >= 1 - z  ] 
#              [ 2x + y >= 1 - z  ] 
#              [ epsilon >= w     ]
#
class MyIneq(Optizelle.VectorValuedFunction):

    # Read in the amount of allowable infeasibility into the problem
    def __init__(self,epsilon_):
        self.epsilon=epsilon_

    # z=h(x) 
    def eval(self,x,z):
        z[0]=x[0]+2.*x[1]+x[2]-1.
        z[1]=2.*x[0]+x[1]+x[2]-1.
        z[2]=self.epsilon-x[3]

    # z=h'(x)dx
    def p(self,x,dx,z):
        z[0]= dx[0]+2.*dx[1]+dx[2]
        z[1]= 2.*dx[0]+dx[1]+dx[2]
        z[2]= -dx[3]

    # xhat=h'(x)*dz
    def ps(self,x,dz,xhat):
        xhat[0]= dz[0]+2.*dz[1]
        xhat[1]= 2.*dz[0]+dz[1]
        xhat[2]= dz[0]+dz[1]
        xhat[3]= -dz[2]

    # xhat=(h''(x)dx)*dz
    def pps(self,x,dx,dz,xhat):
        xhat.fill(0.)

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.stderr.write("Usage: python simple_infeasible_inequality.py <parameters>.json\n")
    raise ValueError("Parameters JSON file required.")

fname=sys.argv[1]

# Set the amount of infeasibility that we want to allow
epsilon = 1e-8

# Generate an initial guess for the primal
x = numpy.array([0.,0.,5.,-5.])

# Generate a vector for the equality multiplier 
y = numpy.array([0.])

# Generate a vector for the inequality multiplier
z = numpy.array([0.,0.,0.])

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
fns.h=MyIneq(epsilon)

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
