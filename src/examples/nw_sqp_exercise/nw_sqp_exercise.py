# Exercise 18.2 in Numerical Optimization by Nocedal and Wright.  This has
# an optimal solution of x = (-1.71,1.59,1.82.-0.763,-0.763).

import Optizelle 
import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions
import Optizelle.EqualityConstrained.Algorithms
import Optizelle.json.EqualityConstrained
import numpy
import sys
from math import exp

# Squares its input
sq = lambda x:x*x

# Indexing for vectors
def itok(i):
    return i-1;

# Indexing for packed storage
def ijtokp(i,j):
    if i>j:
        tmp=i
        i=j
        j=tmp
    return (i-1)+j*(j-1)/2
    
# Indexing function for dense matrices 
def ijtok(i,j,m):
    return (i-1)+(j-1)*m

# Indexing function for dense tensors 
def ijktol(i,j,k,m,n):
    return (i-1)+(j-1)*m+(k-1)*m*n

#
# f(x) = exp(x1 x2 x3 x4 x5) - (1/2) (x1^3 + x2^3 + 1)^2
#
class MyObj(Optizelle.ScalarValuedFunction):

    # Evaluation 
    def eval(self,x):
        return (exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - sq(pow(x[itok(1)],3)+pow(x[itok(2)],3)+1.)/2.)

    # Gradient
    def grad(self,x,g):
        g[itok(1)]= (x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - 3.*sq(x[itok(1)])
            *(pow(x[itok(1)],3) + pow(x[itok(2)],3)+1.))
        g[itok(2)]= (x[itok(1)]*x[itok(3)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
            - 3.*sq(x[itok(2)])
            * (pow(x[itok(1)],3) + pow(x[itok(2)],3) + 1.))
        g[itok(3)]= (x[itok(1)]*x[itok(2)]*x[itok(4)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))
        g[itok(4)]= (x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(5)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))
        g[itok(5)] = (x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))

    # Hessian-vector product
    def hessvec(self,x,dx,H_dx):
        # Allocate memory for the dense Hessian in packed storage
        H=numpy.empty(15)

        # Compute the dense Hessian
        H[ijtokp(1,1)] = (
            sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                - 6.*x[itok(1)]*(pow(x[itok(1)],3)
                + pow(x[itok(2)],3)+1.)
                - 9.*pow(x[itok(1)],4))
        H[ijtokp(1,2)]= (
            x[itok(1)]*x[itok(2)]*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(3)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)])
                - 9.*sq(x[itok(1)])*sq(x[itok(2)]))
        H[ijtokp(1,3)]= (
            x[itok(1)]*sq(x[itok(2)])*x[itok(3)]*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(1,4)]= (
            x[itok(1)]*sq(x[itok(2)])*sq(x[itok(3)])*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(3)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(1,5)]= (
            x[itok(1)]*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(2)]*x[itok(3)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(2,2)] = (
            sq(x[itok(1)])*sq(x[itok(3)])*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                - 6.*x[itok(2)]
                * (pow(x[itok(1)],3)+pow(x[itok(2)],3)+1.)
                - 9.*pow(x[itok(2)],4))
        H[ijtokp(2,3)]= (
            sq(x[itok(1)])*x[itok(2)]*x[itok(3)]*sq(x[itok(4)])*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(4)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(2,4)]= (
            sq(x[itok(1)])*x[itok(2)]*sq(x[itok(3)])*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(3)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(2,5)]= (
            sq(x[itok(1)])*x[itok(2)]*sq(x[itok(3)])*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(3)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(3,3)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(4)])*sq(x[itok(5)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(3,4)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*x[itok(3)]*x[itok(4)]*sq(x[itok(5)])
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(5)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(3,5)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*x[itok(3)]*sq(x[itok(4)])*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(4)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(4,4)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(5)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(4,5)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*x[itok(4)]*x[itok(5)]
                * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)])
                + x[itok(1)]*x[itok(2)]*x[itok(3)]*exp(x[itok(1)]*x[itok(2)]
                * x[itok(3)]*x[itok(4)]*x[itok(5)]))
        H[ijtokp(5,5)]= (
            sq(x[itok(1)])*sq(x[itok(2)])*sq(x[itok(3)])*sq(x[itok(4)])
            * exp(x[itok(1)]*x[itok(2)]*x[itok(3)]*x[itok(4)]*x[itok(5)]))

        # Compute the Hessian-vector product
        H_dx.fill(0.)
        for i in xrange(1,6):
            for j in xrange(1,6):
                H_dx[i-1] += H[ijtokp(i,j)]*dx[j-1]

#
# g(x)= [ x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 ]
#       [ x2 x3 - 5 x4 x5                       ]
#       [ x1^3 + x2^3 + 1                       ]
#
class MyEq(Optizelle.VectorValuedFunction):

    # y=g(x) 
    def eval(self,x,y): 
        y[itok(1)] = (sq(x[itok(1)]) + sq(x[itok(2)]) + sq(x[itok(3)])
            + sq(x[itok(4)]) + sq(x[itok(5)]) - 10.)
        y[itok(2)] = (x[itok(2)]*x[itok(3)] - 5.*x[itok(4)]*x[itok(5)])
        y[itok(3)] = (pow(x[itok(1)],3) + pow(x[itok(2)],3) + 1.)

    # Generate a dense version of the Jacobian
    def generateJac(self,x,jac):
        jac[ijtok(1,1,3)] = 2.*x[itok(1)]
        jac[ijtok(1,2,3)] = 2.*x[itok(2)]
        jac[ijtok(1,3,3)] = 2.*x[itok(3)]
        jac[ijtok(1,4,3)] = 2.*x[itok(4)]
        jac[ijtok(1,5,3)] = 2.*x[itok(5)]
        
        jac[ijtok(2,2,3)] = x[itok(3)]
        jac[ijtok(2,3,3)] = x[itok(2)]
        jac[ijtok(2,4,3)] = -5.*x[itok(5)]
        jac[ijtok(2,5,3)] = -5.*x[itok(4)]
        
        jac[ijtok(3,1,3)] = 3.*sq(x[itok(1)])
        jac[ijtok(3,2,3)] = 3.*sq(x[itok(2)])

    # y=g'(x)dx
    def p(self,x,dx,y):
        # Generate a dense matrix that holds the Jacobian
        jac = numpy.zeros(15)

        # Compute a dense form of the Jacobian
        self.generateJac(x,jac);

        # Compute the Jacobian-vector product
        y.fill(0.)
        for i in xrange(1,4):
            for j in xrange(1,6):
                y[itok(i)] += jac[ijtok(i,j,3)]*dx[itok(j)]

    # z=g'(x)*dy
    def ps(self,x,dy,z):
        # Generate a dense matrix that holds the Jacobian
        jac = numpy.zeros(15)

        # Compute a dense form of the Jacobian
        self.generateJac(x,jac);

        # Compute the Jacobian transpose-vector product
        z.fill(0.)
        for i in xrange(1,4):
            for j in xrange(1,6):
                z[itok(j)] += jac[ijtok(i,j,3)]*dy[itok(i)]

    # z=(g''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        # Generate a dense tensor that holds the second derivative adjoint
        D = numpy.zeros(75)
        D[ijktol(1,1,1,3,5)] = 2.
        D[ijktol(1,2,2,3,5)] = 2.
        D[ijktol(1,3,3,3,5)] = 2.
        D[ijktol(1,4,4,3,5)] = 2.
        D[ijktol(1,5,5,3,5)] = 2.
        
        D[ijktol(2,2,3,3,5)] = 1.
        D[ijktol(2,3,2,3,5)] = 1.
        D[ijktol(2,4,5,3,5)] = -5.
        D[ijktol(2,5,4,3,5)] = -5.

        D[ijktol(3,1,1,3,5)] = 6.*x[itok(1)]
        D[ijktol(3,2,2,3,5)] = 6.*x[itok(2)]

        # Compute the action of this operator on our directions
        z.fill(0.)
        for i in xrange(1,4):
            for j in xrange(1,6):
                for k in xrange(1,6):
                    z[itok(k)] += D[ijktol(i,j,k,3,5)]*dx[itok(j)]*dy[itok(i)]


# Read in the name for the input file
if len(sys.argv)!=2:
    sys.exit("python nw_sqp_exercise.py <parameters>")
fname = sys.argv[1]

# Generate an initial guess for the primal
x = numpy.array([-1.8,1.7,1.9,-0.8,-0.8])

# Generate an initial guess for the dual
y = numpy.zeros(3)

# Create an optimization state
state=Optizelle.EqualityConstrained.State.t(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),x,y)

# Read the parameters from file
Optizelle.json.EqualityConstrained.read(Optizelle.Rm,Optizelle.Rm,
    Optizelle.Messaging(),fname,state)

# Create the bundle of functions 
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
print "The optimal point is:"
for i in xrange(1,6):
    if i==1:
        sys.stdout.write("[ ")
    else:
        sys.stdout.write("  ")
    sys.stdout.write("%13e" % state.x[itok(i)])
    if i==5:
        print " ]"
    else:
        print " ;"

# Write out the final answer to file
Optizelle.json.EqualityConstrained.write_restart(
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging(),"solution.json",state)
