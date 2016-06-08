# In this example, we setup and minimize the Rosenbrock function.

import Optizelle 
import numpy
import sys

#---Objective0---
# Squares its input
sq = lambda x:x*x

# Define the Rosenbrock function where
# 
# f(x,y)=(1-x)^2+100(y-x^2)^2
#
class Rosenbrock(Optizelle.ScalarValuedFunction):
    # Evaluation of the Rosenbrock function
    def eval(self,x):
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]))

    # Gradient
    def grad(self,x,grad): 
        grad[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0])
        grad[1]=200*(x[1]-sq(x[0]))

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
    	H_dx[0] = (1200*sq(x[0])-400*x[1]+2)*dx[0]-400*x[0]*dx[1]
        H_dx[1] = -400*x[0]*dx[0] + 200*dx[1]
#---Objective1---

#---Preconditioner0---
# Define a perfect preconditioner for the Hessian
class RosenHInv(Optizelle.Operator):
    def eval(self,state,dx,result):
        x = state.x
        one_over_det=1./(80000.*sq(x[0])-80000.*x[1]+400.)
        result[0]=one_over_det*(200.*dx[0]+400.*x[0]*dx[1])
        result[1]=(one_over_det*
            (400.*x[0]*dx[0]+(1200.*x[0]*x[0]-400.*x[1]+2.)*dx[1]))
#---Preconditioner1---
    
# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("python rosenbrock.py <parameters>")
fname = sys.argv[1]

#---State0---
# Generate an initial guess for Rosenbrock
x = numpy.array([-1.2,1.0])

# Create an unconstrained state based on this vector
state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x)
#---State1---
    
#---Parameters0---
# Read the parameters from file
Optizelle.json.Unconstrained.read(Optizelle.Rm,fname,state)
#---Parameters1---

#---Functions0---
# Create the bundle of functions 
fns=Optizelle.Unconstrained.Functions.t()
fns.f=Rosenbrock()
fns.PH=RosenHInv()
#---Functions1---

#---Solver0---
# Solve the optimization problem
Optizelle.Unconstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Messaging.stdout,fns,state)
#---Solver1---

#---Extract0---
# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])
#---Extract1---

# Write out the final answer to file
Optizelle.json.Unconstrained.write_restart(Optizelle.Rm,"solution.json",state)
