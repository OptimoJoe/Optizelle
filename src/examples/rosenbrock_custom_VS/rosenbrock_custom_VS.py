# This example minimizes the Rosenbrock function from scratch.  Meaning,
# it runs through a complete example from defining a valid vector space,
# to setting parameters, to solving the problem.  For reference, the optimal
# solution to the Rosenbrock function is (1,1).

import Optizelle 
import Optizelle.Unconstrained.State
import Optizelle.Unconstrained.Functions
import Optizelle.Unconstrained.Algorithms
import Optizelle.json.Unconstrained
import numpy
import sys
import copy

# Defines the vector space used for optimization.
class MyVS(object):
    @staticmethod
    def init(x):
        """Memory allocation and size setting"""
        return copy.deepcopy(x)

    @staticmethod
    def copy(x,y):
        """y <- x (Shallow.  No memory allocation.)"""
        y[:]=x[:]

    @staticmethod
    def scal(alpha,x):
        """x <- alpha * x"""
        for i in xrange(0,len(x)):
            x[i]=alpha*x[i]

    @staticmethod
    def zero(x):
        """x <- 0"""
        for i in xrange(0,len(x)):
            x[i]=0.

    @staticmethod
    def axpy(alpha,x,y):
        """y <- alpha * x + y"""
        for i in xrange(0,len(x)):
            y[i]=alpha*x[i]+y[i]

    @staticmethod
    def innr(x,y):
        """<- <x,y>"""
        return reduce(lambda z,xy:xy[0]*xy[1]+z,zip(x,y),0.)

    @staticmethod
    def rand(x):
        """x <- random"""
        for i in xrange(0,len(x)):
            x[i]=random.uniform(0.,1.) 

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

# Generate an initial guess for Rosenbrock
x = numpy.array([-1.2,1.0])

# Create an unconstrained state based on this vector
state=Optizelle.Unconstrained.State.t(MyVS,Optizelle.Messaging(),x)

# Setup some algorithmic parameters

if 1:
    # Trust-Region Newton's method
    state.H_type = Optizelle.Operators.UserDefined
    state.iter_max = 50;
    state.eps_krylov = 1e-10;

if 0:
    # BFGS
    state.algorithm_class = Optizelle.AlgorithmClass.LineSearch
    state.dir = Optizelle.LineSearchDirection.BFGS
    state.stored_history = 10;
    state.iter_max = 100;

if 0:
    # Newton-CG 
    state.algorithm_class = Optizelle.AlgorithmClass.LineSearch
    state.dir = Optizelle.LineSearchDirection.NewtonCG
    state.H_type = Optizelle.Operators.UserDefined
    state.eps_krylov = 1e-2
    state.iter_max = 50

# Create the bundle of functions 
fns=Optizelle.Unconstrained.Functions.t()
fns.f=Rosenbrock()

# Solve the optimization problem
Optizelle.Unconstrained.Algorithms.getMin(MyVS,Optizelle.Messaging(),fns,state)

# Print out the reason for convergence
print "The algorithm converged due to: %s" % (
    Optizelle.StoppingCondition.to_string(state.opt_stop))

# Print out the final answer
print "The optimal point is: (%e,%e)" % (state.x[0],state.x[1])
