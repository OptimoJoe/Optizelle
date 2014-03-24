# This example minimizes the Rosenbrock function from scratch.  Meaning,
# it runs through a complete example from defining a valid vector space,
# to setting parameters, to solving the problem.  For reference, the optimal
# solution to the Rosenbrock function is (1,1).

from __future__ import print_function
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

    @staticmethod
    def prod(x,y,z):
        """Jordan product, z <- x o y"""
        numpy.copyto(z,x*y)

    @staticmethod
    def id(x):
        """Identity element, x <- e such that x o e = x"""
        x.fill(1.)

    @staticmethod
    def linv(x,y,z):
        """Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y"""
        numpy.copyto(z,numpy.divide(y,x))

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

# Define a perfect preconditioner for the Hessian
class RosenHInv(Optizelle.Operator):
    def eval(self,state,dx,result):
        x = state.x
        one_over_det=1./(80000.*sq(x[0])-80000.*x[1]+400.)
        result[0]=one_over_det*(200.*dx[0]+400.*x[0]*dx[1])
        result[1]=(one_over_det*
            (400.*x[0]*dx[0]+(1200.*x[0]*x[0]-400.*x[1]+2.)*dx[1]))

# Define a custom messaging object
class MyMessaging(Optizelle.Messaging):
    """Defines how we output messages to the user"""
   
    def print(self,msg):
        """Prints out normal diagnostic information"""
        sys.stdout.write("PRINT:  %s\n" %(msg))

    def error(self,msg):
        """Prints out error information"""
        sys.stderr.write("ERROR:  %s\n" %(msg))

def serialize_MyVS(x):
    """Serializes an array for the vector space MyVS""" 

    # Create the json representation
    x_json="[ "
    for i in xrange(x.size):
        x_json  += str(x[i]) + ", "
    x_json=x_json[0:-2]
    x_json +=" ]"

    return x_json

def deserialize_MyVS(x,x_json):
    """Deserializes an array for the vector space MyVS""" 

    # Eliminate all whitespace
    x_json="".join(x_json.split())

    # Check if we're a vector
    if x_json[0:1]!="[" or x_json[-1:]!="]":
        raise TypeError("Attempted to deserialize a non-numpy.array vector.")

    # Eliminate the initial and final delimiters
    x_json=x_json[1:-1]

    # Create a list of the numbers involved 
    x_json=x_json.split(",")

    # Convert the strings to numbers
    x_json=map(lambda x:float(x),x_json)

    # Create a MyVS vector
    return array.array('d',x_json)

# Define a state manipulator that writes out the optimization state at
# each iteration.
class MyRestartManipulator(Optizelle.StateManipulator):
    def eval(self,fns,state,loc):
        # At the end of the optimization iteration, write the restart file
        if loc == Optizelle.OptimizationLocation.EndOfOptimizationIteration:
            # Create a reasonable file name
            ss = "rosenbrock_advanced_api_%04d.json" % (state.iter)
                
            # Write the restart file
            Optizelle.json.Unconstrained.write_restart( 
               MyVS,MyMessaging(),ss,state)
    
# Read in the name for the input file
if not(len(sys.argv)==2 or len(sys.argv)==3):
    sys.exit("python rosenbrock_advanced_api.py <parameters>\n" +
             "python rosenbrock_advanced_api.py <parameters> <restart>")
pname = sys.argv[1]
rname = sys.argv[2] if len(sys.argv)==3 else ""

# Generate an initial guess for Rosenbrock
x = numpy.array([-1.2,1.0])

# Create an unconstrained state based on this vector
state=Optizelle.Unconstrained.State.t(MyVS,MyMessaging(),x)

# If we have a restart file, read in the parameters 
if len(sys.argv)==3:
    Optizelle.json.Unconstrained.read_restart(
        MyVS,MyMessaging(),rname,x,state)

# Read the parameters from file
Optizelle.json.Unconstrained.read(Optizelle.Rm,Optizelle.Messaging(),
    pname,state)

# Create the bundle of functions 
fns=Optizelle.Unconstrained.Functions.t()
fns.f=Rosenbrock()
fns.PH=RosenHInv()

# Solve the optimization problem
Optizelle.Unconstrained.Algorithms.getMin(
    MyVS,MyMessaging(),MyRestartManipulator(),fns,state)

# Print out the reason for convergence
print("The algorithm converged due to: %s" % (
    Optizelle.StoppingCondition.to_string(state.opt_stop)))

# Print out the final answer
print("The optimal point is: (%e,%e)" % (state.x[0],state.x[1]))

# Write out the final answer to file
Optizelle.json.Unconstrained.write_restart(
    Optizelle.Rm,Optizelle.Messaging(),"solution.json",state)
