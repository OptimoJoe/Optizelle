# In this example, we setup and minimize the Rosenbrock function

# Load in peopt
import peopt
from math import *

# Create a vector space 
class ListVS:
    def copy(self,x):
        return x
    def scal(self,alpha,x):
        return map(lambda y:alpha*y,x)
    def zero(self,x):
        return [0]*len(x)
    def axpy(self,alpha,x,y):
        return map(lambda z:alpha*z[0]+z[1],zip(x,y))
    def innr(self,x,y):
        return reduce(lambda x,y:x+y[0]*y[1],zip(x,y),0)

# Define the Rosenbrock function where
# 
# f(x,y)=(1-x)^2+100(y-x^2)^2
#
sq = lambda x:x*x
class Rosenbrock:
    # Evaluation of the Rosenbrock function
    def eval(self,x):
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]))

    # Gradient
    def grad(self,x): 
        g=[0]*len(x)
        g[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0])
        g[1]=200*(x[1]-sq(x[0]))
        return g

    # Hessian-vector product
    def hessvec(self,x,dx):
        H_dx = [0]*len(x)
    	H_dx[0] = (1200*sq(x[0])-400*x[1]+2)*dx[0]-400*x[0]*dx[1]
        H_dx[1] = -400*x[0]*dx[0] + 200*dx[1]
        return H_dx

# Create a bundle of vector spaces 
class vs:
    X = ListVS()

# Create a bundle of functions
class fns:
    f = Rosenbrock()

# Create a bundle of points
class pts:
    x = [-1.2,1]

x_sol = peopt.getMin(vs(),fns(),pts(),'run.peopt')
#print x_sol
