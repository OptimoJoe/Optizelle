# Optimize a simple optimization problem with an optimal solution 
# of (2-sqrt(2)/2,2-sqrt(2)/2).

import Optizelle 
import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions
import Optizelle.EqualityConstrained.Algorithms
import Optizelle.json.EqualityConstrained
import numpy
import sys

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

    # z=g'(x)*dy
    def ps(self,x,dy,z):
        z[0] = 2.*(x[0]-2.)*dy[0]
        z[1] = 2.*(x[1]-2.)*dy[0]

    # z=(g''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        z[0] = 2.*dx[0]*dy[0];
        z[1] = 2.*dx[1]*dy[0];
#---EqualityConstraint1---

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("simple_equality.py <parameters>")
fname=sys.argv[1]

# Generate an initial guess 
x = numpy.array([2.1,1.1])

# Allocate memory for the equality multiplier 
y = numpy.array([0.])

# Create an optimization state
state=Optizelle.EqualityConstrained.State.t(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),x,y)

# Read the parameters from file
Optizelle.json.EqualityConstrained.read(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),fname,state)

# Create a bundle of functions
fns=Optizelle.EqualityConstrained.Functions.t()
fns.f=MyObj()
fns.g=MyEq()

# Solve the optimization problem
Optizelle.EqualityConstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.StoppingCondition.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])

# Write out the final answer to file
Optizelle.json.EqualityConstrained.write_restart(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),"solution.json",state)
