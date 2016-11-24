# This example minimizes a simple quadratic function.  For reference, the 
# optimal solution to this function is (1,1,1).

import Optizelle 
import numpy
import sys

# Squares its input
sq = lambda x:x*x

# Define the quadratic function 
# 
# f(x,y,z)=(x-1)^2+(2y-2)^2+(3z-3)^2
#
class Quad(Optizelle.ScalarValuedFunction):

    # Evaluation of the quadratic function
    def eval(self,x):
        return sq(x[0]-1.)+sq(2*x[1]-2.)+sq(3*x[2]-3.)

    # Gradient
    def grad(self,x,g):
        g[0]=2*x[0]-2
        g[1]=8*x[1]-8
        g[2]=18*x[2]-18

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
    	H_dx[0]= 2*dx[0]
        H_dx[1]= 8*dx[1] 
        H_dx[2]= 18*dx[2] 

# Define an almost perfect preconditioner for the Hessian
class QuadHInv(Optizelle.Operator):
    def eval(self,state,dx,result):
        result[0]=dx[0]/2.
        result[1]=dx[1]/8.
        result[2]=dx[2]

# Read in the name for the input file
if len(sys.argv)!=2:
    sys.stderr.write("Usage: python quadratic.py <parameters>.json\n")
    raise ValueError("Parameters JSON file required.")

fname = sys.argv[1]

# Generate an initial guess 
x = numpy.array([-1.2,1.1,2.])

# Create an unconstrained state based on this vector
state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x)

# Read the parameters from file
Optizelle.json.Unconstrained.read(Optizelle.Rm,fname,state)

# Create the bundle of functions 
fns=Optizelle.Unconstrained.Functions.t()
fns.f=Quad()
fns.PH=QuadHInv()

# Solve the optimization problem
Optizelle.Unconstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e,%e)" % (state.x[0],state.x[1],state.x[2])

# Write out the final answer to file
Optizelle.json.Unconstrained.write_restart(Optizelle.Rm,"solution.json",state)
