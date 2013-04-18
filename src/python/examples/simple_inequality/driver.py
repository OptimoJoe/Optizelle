# In this example, we setup and minimize the problem
# min (x+1)^2+(y+1)^2 st x + 2y >= 1, 2x + y >=1

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
    def prod(self,x,y):
        return map(lambda z:z[0]*z[1],zip(x,y))
    def id(self,x):
        return [1]*len(x)
    def linv(self,x,y):
        return map(lambda z:z[1]/z[0],zip(x,y))
    def barr(self,x):
        return reduce(lambda x,y:x+log(y),x,0)
    def srch(self,x,y):
        alpha = float("inf")
        for i in xrange(0,len(x)):
            if x[i] < 0:
                alpha0 = -y[i]/x[i]
                if alpha0 < alpha:
                    alpha=alpha0
        return alpha
    def symm(self,x):
        return x 

# Define the objective 
# 
# f(x,y)=(x+1)^2+(y+1)^2
#
sq = lambda x:x*x
class QuadObj:
    # Evaluation of the Rosenbrock function
    def eval(self,x):
        return sq(x[0]+1)+sq(x[1]+1) 

    # Gradient
    def grad(self,x): 
        g=[0]*len(x)
        g[0]=2*(x[0]+1)
        g[1]=2*(x[1]+1)
        return g

    # Hessian-vector product
    def hessvec(self,x,dx):
        H_dx = [0]*len(x)
    	H_dx[0] = 2*dx[0] 
        H_dx[1] = 2*dx[1] 
        return H_dx

# Define the constraints 
#
# g(x,y)= [ x + 2y - 1] 
#         [ 2x + y - 1] 
#
class LinConst:
    # y=g(x) 
    def eval(self,x):
        y = [0]*2
        y[0]=x[0]+2*x[1]-1
        y[1]=2*x[0]+x[1]-1
        return y

    # y=g'(x)dx
    def p(self,x,dx):
        y = [0]*2
        y[0]= dx[0] + 2*dx[1]
        y[1]= 2*dx[0] + dx[1] 
        return y

    # z=g'(x)*dy
    def ps(self,x,dy):
        z = [0]*2
        z[0]= dy[0] + 2*dy[1] 
        z[1]= 2*dy[0] + dy[1]
        return z

    # z=(g''(x)dx)*dy
    def pps(self,x,dx,dy):
        z = [0]*2
        z[0] = 0.
        z[1] = 0. 
        return z

# Create a bundle of vector spaces 
class vs:
    X = ListVS()
    Z = ListVS()

# Create a bundle of functions
class fns:
    f = QuadObj()
    h = LinConst()

# Create a bundle of points
class pts:
    x = [3.,4.]
    dx = [5.,6.]
    dxx = [7.,8.]
    z = [1.,1.]
    dz = [9.,10.]

#peopt.diagnostics(vs(),fns(),pts())
sol = peopt.getMin(vs(),fns(),pts(),'simple_inequality.peopt')
print sol[0]
