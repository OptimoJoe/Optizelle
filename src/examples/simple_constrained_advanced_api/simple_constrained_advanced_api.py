# Optimize a simple optimization problem with an optimal solution
# of (1/3,1/3)

from __future__ import print_function
import Optizelle 
import Optizelle.Constrained.State
import Optizelle.Constrained.Functions
import Optizelle.Constrained.Algorithms
import Optizelle.json.Constrained
import Optizelle.json.Serialization
import sys
import copy
import array
import math

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

    @staticmethod
    def prod(x,y,z):
        """Jordan product, z <- x o y"""
        for i in xrange(0,len(x)):
            z[i]=x[i]*y[i]

    @staticmethod
    def id(x):
        """Identity element, x <- e such that x o e = x"""
        for i in xrange(0,len(x)):
            x[i]=1.

    @staticmethod
    def linv(x,y,z):
        """Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y"""
        for i in xrange(0,len(x)):
            z[i]=y[i]/x[i]

    @staticmethod
    def barr(x):
        """Barrier function, <- barr(x) where x o grad barr(x) = e"""
        return reduce(lambda x,y:x+math.log(y),x,0.)
        
    @staticmethod
    def srch(x,y):
        """Line search, <- argmax {alpha \in Real >= 0 : alpha x + y >= 0} where y > 0"""
        alpha = float("inf")
        for i in xrange(0,len(x)):
            if x[i] < 0:
                alpha0 = -y[i]/x[i]
                if alpha0 < alpha:
                    alpha=alpha0
        return alpha

    @staticmethod
    def symm(x):
        """Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric operator"""
        pass

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

    # z=g'(x)*dy
    def ps(self,x,dy,z):
        z[0]= dy[0]
        z[1]= 2.*dy[0]

    # z=(g''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        MyVS.zero(z)

# Define simple inequalities 
#
# h(x,y)= [ 2x + y >= 1 ] 
#
class MyIneq(Optizelle.VectorValuedFunction):

    # y=h(x) 
    def eval(self,x,y):
        y[0]=2.*x[0]+x[1]-1.

    # y=h'(x)dx
    def p(self,x,dx,y):
        y[0]= 2.*dx[0]+dx[1]

    # z=h'(x)*dy
    def ps(self,x,dy,z):
        z[0]= 2.*dy[0]
        z[1]= dy[0]

    # z=(h''(x)dx)*dy
    def pps(self,x,dx,dy,z):
        MyVS.zero(z)

#---Serialization0---
def serialize_MyVS(x,name,iter):
    """Serializes an array for the vector space MyVS""" 

    # Create the filename where we put our vector
    fname = "./restart/%s.%04d.txt" % (name,iter)

    # Actually write the vector there
    fout = open(fname,"w");
    for i in xrange(0,len(x)):
        fout.write("%1.16e\n" % x[i]) 

    # Close out the file
    fout.close()

    # Use this filename as the json string 
    x_json = "\"%s\"" % fname
    return x_json

def deserialize_MyVS(x_,x_json):
    """Deserializes an array for the vector space MyVS""" 

    # Eliminate all whitespace
    x_json="".join(x_json.split())

    # Eliminate the initial and final delimiters
    x_json=x_json[1:-1]

    # Open the file for reading
    fin = open(x_json,"r")

    # Allocate a new vector to return
    x = copy.deepcopy(x_)
    
    # Read in each of the elements
    for i in xrange(0,len(x)):
        x[i] = float(fin.readline())

    # Close out the file
    fin.close()

    # Return the result 
    return x 

# Register the serialization routines for arrays 
def MySerialization():
    Optizelle.json.Serialization.serialize.register(
        serialize_MyVS,array.array)
    Optizelle.json.Serialization.deserialize.register(
        deserialize_MyVS,array.array)
#---Serialization1---

# Define a state manipulator that writes out the optimization state at
# each iteration.
class MyRestartManipulator(Optizelle.StateManipulator):
    def eval(self,fns,state,loc):
        # At the end of the optimization iteration, write the restart file
        if loc == Optizelle.OptimizationLocation.EndOfOptimizationIteration:
            # Create a reasonable file name
            ss = "simple_constrained_advanced_api_%04d.json" % (state.iter)
                
            # Write the restart file
            Optizelle.json.Constrained.write_restart( 
               MyVS,MyVS,MyVS,Optizelle.Messaging(),ss,state)

# Register the serialization routines
MySerialization()

# Read in the name for the input file
if not(len(sys.argv)==2 or len(sys.argv)==3):
    sys.exit("python simple_constrained_advanced_api.py <parameters>\n" +
             "python simple_constrained_advanced_api.py <parameters> <restart>")
pname = sys.argv[1]
rname = sys.argv[2] if len(sys.argv)==3 else ""

# Generate an initial guess 
x = array.array('d',[2.1,1.1])

# Allocate memory for the equality multiplier 
y = array.array('d',[0.])

# Allocate memory for the inequality multiplier 
z = array.array('d',[0.])

# Create an optimization state
state=Optizelle.Constrained.State.t(
    MyVS,MyVS,MyVS,Optizelle.Messaging(),x,y,z)

# If we have a restart file, read in the parameters 
if len(sys.argv)==3:
    Optizelle.json.Constrained.read_restart(
        MyVS,MyVS,MyVS,Optizelle.Messaging(),rname,x,y,z,state)

# Read the parameters from file
Optizelle.json.Constrained.read(
    MyVS,MyVS,MyVS,Optizelle.Messaging(),pname,state)

# Create a bundle of functions
fns=Optizelle.Constrained.Functions.t()
fns.f=MyObj()
fns.g=MyEq()
fns.h=MyIneq()

# Solve the optimization problem
Optizelle.Constrained.Algorithms.getMin(
    MyVS,MyVS,MyVS,Optizelle.Messaging(),fns,state,MyRestartManipulator())

# Print out the reason for convergence
print("The algorithm converged due to: %s" % (
    Optizelle.OptimizationStop.to_string(state.opt_stop)))

# Print out the final answer
print("The optimal point is: (%e,%e)" % (state.x[0],state.x[1]))

# Write out the final answer to file
Optizelle.json.Constrained.write_restart(MyVS,MyVS,MyVS,
    Optizelle.Messaging(),"solution.json",state)
