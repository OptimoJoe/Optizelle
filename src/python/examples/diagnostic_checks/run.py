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


# Define some utility function where
#
# g(x,y)= [ cos(x)sin(y) ]
#         [ 3 x^2 y + y^3]
#         [ log(x) + 3y^5]
#
cub = lambda x:x*x*x
quad = lambda x:x*x*x*x
quint = lambda x:x*x*x*x*x
class Utility:
    # y=f(x) 
    def eval(self,x):
        y = [0]*3
        y[0]=cos(x[0])*sin(x[1])
        y[1]=3.*sq(x[0])*x[1]+cub(x[1])
        y[2]=log(x[0])+3.*quint(x[1])
        return y

    # y=f'(x)dx
    def p(self,x,dx):
        y = [0]*3
        y[0]= -sin(x[0])*sin(x[1])*dx[0] \
              +cos(x[0])*cos(x[1])*dx[1]
        y[1]= 6.*x[0]*x[1]*dx[0] \
              +(3.*sq(x[0])+3.*sq(x[1]))*dx[1]
        y[2]= 1./x[0]*dx[0] \
              +15.*quad(x[1])*dx[1]
        return y

    # z=f'(x)*dy
    def ps(self,x,dy):
        z = [0]*2
        z[0]= -sin(x[0])*sin(x[1])*dy[0] \
              +6.*x[0]*x[1]*dy[1] \
              +1./x[0]*dy[2]
        z[1]= cos(x[0])*cos(x[1])*dy[0] \
              +(3.*sq(x[0])+3.*sq(x[1]))*dy[1] \
              +15.*quad(x[1])*dy[2]
        return z

    # z=(f''(x)dx)*dy
    def pps(self,x,dx,dy):
        z = [0]*2
        z[0] = (-cos(x[0])*dx[0]*sin(x[1])-sin(x[0])*cos(x[1])*dx[1])*dy[0] \
               +(6.*dx[0]*x[1] + 6.*x[0]*dx[1])*dy[1] \
               +(-1./sq(x[0])*dx[0])*dy[2]
        z[1] = (-sin(x[0])*dx[0]*cos(x[1])-cos(x[0])*sin(x[1])*dx[1])*dy[0] \
               +(6.*x[0]*dx[0]+6.*x[1]*dx[1])*dy[1] \
               +(60.*cub(x[1])*dx[1])*dy[2]
        return z

# Create a bundle of vector spaces 
class vs:
    X = ListVS()
    Y = ListVS()

# Create a bundle of functions
class fns:
    f = Rosenbrock()
    g = Utility()

# Create a bundle of points
class pts:
    x = [1,2]
    dx = [3/1000.,4/1000.]
    dxx = [5/1000.,6/1000.]
    y = [7,8,9]
    dy = [9,10,11]

peopt.fd_check(vs(),fns(),pts())
