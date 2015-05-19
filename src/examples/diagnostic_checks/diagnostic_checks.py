# This example demonstrates how to run a series of diagnostic tests
# on functions and then immediately exit.

import Optizelle 
import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions
import Optizelle.InequalityConstrained.Algorithms
import Optizelle.json.InequalityConstrained
import numpy
import sys
import math
from math import cos
from math import sin

# Squares its input
sq = lambda x:x*x

# Cubes its input
cub = lambda x:x*x*x

# Quads its input
quad = lambda x:x*x*x*x

# Quints its input
quint = lambda x:x*x*x*x*x

# Defines a log function that returns NaN if its argument is negative.
# Python's built-in log function raises an exception instead, which differs
# from C++ and messes up this example.  Optizelle will correctly compute in the
# presence of NaNs, but it will exit if it detects Python exceptions.
def log(x):
    if x<=0.:
        return float("nan")
    else:
        return math.log(x)

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

# Define some utility function where
#
# g(x)= [ cos(x1) sin(x2)   ]
#       [ 3 x1^2 x2 + x2 ^3 ]
#       [ log(x1) + 3 x2 ^5 ]
#
class Utility(Optizelle.VectorValuedFunction):
    # y=g(x) 
    def eval(self,x,y):
        y[0]=cos(x[0])*sin(x[1])
        y[1]=3.*sq(x[0])*x[1]+cub(x[1])
        y[2]=log(x[0])+3.*quint(x[1])

    # y=g'(x)dx
    def p(self,x,dx,y):
        y[0]= (-sin(x[0])*sin(x[1])*dx[0]
              +cos(x[0])*cos(x[1])*dx[1])
        y[1]= (6.*x[0]*x[1]*dx[0]
              +(3.*sq(x[0])+3.*sq(x[1]))*dx[1])
        y[2]= (1./x[0]*dx[0]
              +15.*quad(x[1])*dx[1])

    # z=g'(x)*dy
    def ps(self,x,dy,z):
        z[0]= (-sin(x[0])*sin(x[1])*dy[0]
              +6.*x[0]*x[1]*dy[1]
              +1./x[0]*dy[2])
        z[1]= (cos(x[0])*cos(x[1])*dy[0]
              +(3.*sq(x[0])+3.*sq(x[1]))*dy[1]
              +15.*quad(x[1])*dy[2])

    # z=(g''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        z[0] = ((-cos(x[0])*dx[0]*sin(x[1])-sin(x[0])*cos(x[1])*dx[1])*dy[0]
               +(6.*dx[0]*x[1] + 6.*x[0]*dx[1])*dy[1]
               +(-1./sq(x[0])*dx[0])*dy[2])
        z[1] = ((-sin(x[0])*dx[0]*cos(x[1])-cos(x[0])*sin(x[1])*dx[1])*dy[0]
               +(6.*x[0]*dx[0]+6.*x[1]*dx[1])*dy[1]
               +(60.*cub(x[1])*dx[1])*dy[2])

# Allocate memory for an initial guess and equality multiplier 
x = numpy.array([1.2,2.3])
z = numpy.zeros(3)

# Create an optimization state
state=Optizelle.InequalityConstrained.State.t(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),x,z)

# Modify the state so that we just run our diagnostics and exit
state.dscheme = Optizelle.DiagnosticScheme.DiagnosticsOnly
state.f_diag = Optizelle.FunctionDiagnostics.SecondOrder
state.x_diag = Optizelle.VectorSpaceDiagnostics.Basic
state.h_diag = Optizelle.FunctionDiagnostics.SecondOrder
state.z_diag = Optizelle.VectorSpaceDiagnostics.EuclideanJordan
state.L_diag = Optizelle.FunctionDiagnostics.SecondOrder

# Create a bundle of functions
fns=Optizelle.InequalityConstrained.Functions.t()
fns.f=Rosenbrock()
fns.h=Utility()

# Even though this looks like we're solving an optimization problem,
# we're actually just going to run our diagnostics and then exit.
Optizelle.InequalityConstrained.Algorithms.getMin(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),fns,state)
