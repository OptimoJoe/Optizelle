# Optimize a simple optimization problem with an optimal solution
# of (1/3,1/3)

import Optizelle 
import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions
import Optizelle.InequalityConstrained.Algorithms
import Optizelle.json.InequalityConstrained
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

# Define simple inequalities 
#
# h(x,y)= [ x + 2y >= 1 ] 
#         [ 2x + y >= 1 ] 
#
class MyIneq(Optizelle.VectorValuedFunction):

    # y=h(x) 
    def eval(self,x,y):
        y[0]=x[0]+2.*x[1]-1.
        y[1]=2.*x[0]+x[1]-1.

    # y=h'(x)dx
    def p(self,x,dx,y):
        y[0]= dx[0]+2.*dx[1]
        y[1]= 2.*dx[0]+dx[1]

    # z=h'(x)*dy
    def ps(self,x,dy,z):
        z[0]= dy[0]+2.*dy[1]
        z[1]= 2.*dy[0]+dy[1]

    # z=(h''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        z.fill(0.)

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("simple_inequality.py <parameters>")
fname=sys.argv[1]

# Generate an initial guess 
x = numpy.array([2.1,1.1])

# Allocate memory for the inequality multiplier 
z = numpy.array([0.,0.])

# Create an optimization state
state=Optizelle.InequalityConstrained.State.t(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),x,z)

# Read the parameters from file
Optizelle.json.InequalityConstrained.read(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),fname,state)

# Create a bundle of functions
fns=Optizelle.InequalityConstrained.Functions.t()
fns.f=MyObj()
fns.h=MyIneq()

# Solve the optimization problem
Optizelle.InequalityConstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.StoppingCondition.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])

# Write out the final answer to file
Optizelle.json.InequalityConstrained.write_restart(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),"solution.json",state)
