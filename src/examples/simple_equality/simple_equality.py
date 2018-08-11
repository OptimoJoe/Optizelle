# Optimize a simple optimization problem with an optimal solution
# of (2-sqrt(2)/2,2-sqrt(2)/2).

import Optizelle
import numpy
import sys

#---Objective0---
# Squares its input
sq = lambda x:x*x

# Define a simple objective where
#
# f(x,y)=x^2+y^2
#
class MyObj(Optizelle.ScalarValuedFunction):

    # Evaluation
    def eval(self,x):
        return sq(x[0])+sq(x[1])

    # Gradient
    def grad(self,x,grad):
        grad[0]=2.*x[0]
        grad[1]=2.*x[1]

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        H_dx[0]=2.*dx[0]
        H_dx[1]=2.*dx[1]
#---Objective1---

#---EqualityConstraint0---
# Define a simple equality constraint
#
# g(x,y)= [ (x-2)^2 + (y-2)^2 = 1 ]
#
class MyEq(Optizelle.VectorValuedFunction):

    # y=g(x)
    def eval(self,x,y):
        y[0] = sq(x[0]-2.)+sq(x[1]-2.)-1.

    # y=g'(x)dx
    def p(self,x,dx,y):
        y[0] = 2.*(x[0]-2.)*dx[0]+2.*(x[1]-2.)*dx[1]

    # xhat=g'(x)*dy
    def ps(self,x,dy,xhat):
        xhat[0] = 2.*(x[0]-2.)*dy[0]
        xhat[1] = 2.*(x[1]-2.)*dy[0]

    # xhat=(g''(x)dx)*dy
    def pps(self,x,dx,dy,xhat):
        xhat[0] = 2.*dx[0]*dy[0]
        xhat[1] = 2.*dx[1]*dy[0]
#---EqualityConstraint1---

#---Preconditioner0---
# Define a Schur preconditioner for the equality constraints
class MyPrecon(Optizelle.Operator):
    def eval(self,state,dy,result):
        result[0]=dy[0]/sq(4.*(x[0]-2.)+4.*sq(x[1]-2.))
#---Preconditioner1---

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("simple_equality.py <parameters>")
fname=sys.argv[1]

#---State0---
# Generate an initial guess
x = numpy.array([2.1,1.1])

# Allocate memory for the equality multiplier
y = numpy.array([0.])

# Create an optimization state
state=Optizelle.EqualityConstrained.State.t(Optizelle.Rm,Optizelle.Rm,x,y)
#---State1---

#---Parameters0---
# Read the parameters from file
Optizelle.json.EqualityConstrained.read(Optizelle.Rm,Optizelle.Rm,fname,state)
#---Parameters1---

#---Functions0---
# Create a bundle of functions
fns=Optizelle.EqualityConstrained.Functions.t()
fns.f=MyObj()
fns.g=MyEq()
fns.PSchur_left=MyPrecon()
#---Functions1---

#---Solver0---
# Solve the optimization problem
Optizelle.EqualityConstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)
#---Solver1---

#---Extract0---
# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])
#---Extract1---

# Write out the final answer to file
Optizelle.json.EqualityConstrained.write_restart(
    Optizelle.Rm,Optizelle.Rm,"solution.json",state)
