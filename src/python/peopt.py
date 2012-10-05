""" The pretty efficient optimization library. """
__author__ = "Joseph Young <josyoun@sandia.gov>"

# Set the exposed functions
__all__ = ['diagnostics','getMin','PeoptError']

# Make sure that we can inspect our functions
import inspect

# Load in the math functions
from math import *

# Make sure we have access to the copy functions
import copy

# Load in our extension
import ctypes
ctypes.cdll.LoadLibrary("libpeopt.so")
libpeopt=ctypes.CDLL("libpeopt.so")

libpeopt.pypeopt.argtypes= \
    [ctypes.py_object,ctypes.py_object,ctypes.py_object,ctypes.py_object,
    ctypes.py_object]
libpeopt.pypeopt.restype = ctypes.py_object 

# Define an exception in the case that we have difficulty with the problem
# setup or running peopt
class PeoptError(Exception):
    """
    The general exception for any errors peopt related.  Its value is
    always a string.
    """
    def __init__(self,value):
        self.value=value
    def __str__(self):
        repr(self.value)

# Validates a variety of objects by making sure they include certain methods
# that require a certain number of arguments
def validate_object(v,req_fns):
    # Check each of the above functions
    for fn in req_fns:
        # First, make sure that each function is a member of the vector space 
        if not hasattr(v,fn[0]):
            raise PeoptError('Missing the ' + fn[0] + ' function.')

        # Next, make sure that each member is actually a function and
        # accepts the correct number of arguments.
        if not hasattr(getattr(v,fn[0]),'__call__') or not\
            len(inspect.getargspec(getattr(v,fn[0])).args) == fn[1]:
            raise PeoptError('Function ' + fn[0] + ' must accept ' + 
                str(fn[1]-1) + ' arguments ('+ str(fn[1]) +
                ' including self.)')

# Insure that the vector space includes the functions copy, scal, zero, axpy,
# and innr.  If not, throw an exception.
def validate_vector_space(v,name):
    req_fns = [('copy',2),('scal',3),('zero',2),('axpy',4),('innr',3)]

    # In the case there's a problem, print out a nice error and then throw
    # the original exception.
    try:               
        validate_object(v,req_fns)  
    except PeoptError as e:
        print('Error in the vector space ' + name + '.\n\n' + e.value + '\n\n' +
            'We require the vector space to be a class of the form:\n' + 
            'class ' + name + ':\n' +
            '    copy(self,x): ...\n' +
            '    scal(self,alpha,x): ...\n' +
            '    zero(self,x): ...\n' +
            '    axpy(self,alpha,x,y): ...\n' +
            '    innr(self,x,y,z): ...\n')
        raise

# Insure that the vector space includes the functions prod, id, linv, srch,
# and barr.  If not, throw an exception
def validate_eja(v,name):
    req_fns = [('copy',2),('scal',3),('zero',2),('axpy',4),('innr',3)]
    req_fns = req_fns+[('prod',3),('id',2),('linv',3),('srch',3),('barr',2)]

    # In the case there's a problem, print out a nice error and then throw
    # the original exception.
    try:               
        validate_object(v,req_fns)  
    except PeoptError as e:
        print('Error in the vector space ' + name + '.\n\n' + e.value + '\n\n' +
            'We require the vector space to be a class of the form:\n' + 
            'class ' + name + ':\n' +
            '    copy(self,x): ...\n' +
            '    scal(self,alpha,x): ...\n' +
            '    zero(self,x): ...\n' +
            '    axpy(self,alpha,x,y): ...\n' +
            '    innr(self,x,y,z): ...\n' +
            '    prod(self,x,y): ...\n' +
            '    id(self,x): ...\n' +
            '    linv(self,x,y): ...\n' +
            '    srch(self,x,y): ...\n' +
            '    barr(self,x): ...\n')
        raise

# Define the different problem types 
class ProblemClass:
    Unconstrained, EqualityConstrained, InequalityConstrained, Constrained\
        = range(4)

# Converts the optimization problem type to a string
def pclass_to_str(pclass):
    if pclass==ProblemClass.Unconstrained:
        return "unconstrained"
    elif pclass==ProblemClass.EqualityConstrained:
        return "equality constrained"
    elif pclass==ProblemClass.InequalityConstrained:
        return "inequality constrained"
    elif pclass==ProblemClass.Constrained:
        return "constrained"

# Determine the type of problem from the vector spaces involved.
def type_from_vs(vs):
    try:
        if not hasattr(vs,'X'):
            raise PeoptError(
                'Every collection of vector spaces must include X.')
        if hasattr(vs,'Y') and hasattr(vs,'Z'):
            return ProblemClass.Constrained
        elif hasattr(vs,'Y'):
            return ProblemClass.EqualityConstrained
        elif hasattr(vs,'Z'):
            return ProblemClass.InequalityConstrained
        else:
            return ProblemClass.Unconstrained
    except PeoptError as e:
            print e.value + '\n'
            raise

# Validate a collection of vector spaces
def validate_vector_spaces(vs):
    # First, get the problem type
    pclass = type_from_vs(vs)

    # First, we always validate X
    validate_vector_space(vs.X,'X')

    # Next, we validate based on the type
    if(pclass == ProblemClass.EqualityConstrained or pclass == ProblemClass.Constrained):
        validate_vector_space(vs.Y,'Y')
    if(pclass == ProblemClass.InequalityConstrained or pclass == ProblemClass.Constrained):
        validate_eja(vs.Z,'Z')

# Validates the objective function.  This means that we need the eval, grad,
# and hessvec functions defined.
def validate_objective(fn,name):
    req_fns = [('eval',2),('grad',2),('hessvec',3)]
    excn = lambda x:raise_(PeoptError(x))

    # In the case there's a problem, print out a nice error and then throw
    # the original exception.
    try:               
        validate_object(fn,req_fns)  
    except PeoptError as e:
        print('Error in the objective function '+ name+'.\n\n' +e.value+'\n\n'+ 
            'We require the objective function to be a class of the form:\n' +
            'class ' + name + ':\n' + 
            '    eval(self,x): ...\n' + 
            '    grad(self,x): ...\n' + 
            '    hessvec(self,x,dx): ...\n')
        raise

# Validates the constraints.  This means that we need the eval, p,
# ps, and pps functions defined.
def validate_constraint(fn,name):
    req_fns = [('eval',2),('p',3),('ps',3),('pps',4)]
    excn = lambda x:raise_(PeoptError(x))

    # In the case there's a problem, print out a nice error and then throw
    # the original exception.
    try:               
        validate_object(fn,req_fns)  
    except PeoptError as e:
        print('Error in the constraint function ' +name+'.\n\n' +e.value+'\n\n'+
            'We require the constraint function to be a class of the form:\n' +
            'class ' + name + ':\n' + 
            '    eval(self,x): ...\n' + 
            '    p(self,x,dx): ...\n' + 
            '    ps(self,x,dy): ...\n' + 
            '    pps(self,x,dx,dy): ...\n') 
        raise

# Determine the type of problem from the functions involved.
def type_from_fns(fns):
    try:
        if not hasattr(fns,'f'):
            raise PeoptError('Every collection of functions must include f.')
        if hasattr(fns,'g') and hasattr(fns,'h'):
            return ProblemClass.Constrained
        elif hasattr(fns,'g'):
            return ProblemClass.EqualityConstrained
        elif hasattr(fns,'h'):
            return ProblemClass.InequalityConstrained
        else:
            return ProblemClass.Unconstrained
    except PeoptError as e:
        print e.value + '\n'
        raise

# Validate a collection of functions 
def validate_functions(fns):
    # First, get the problem type
    pclass = type_from_fns(fns)

    # First, we always validate g
    validate_objective(fns.f,'g')

    # Next, we validate based on the type
    if(pclass == ProblemClass.EqualityConstrained or pclass == ProblemClass.Constrained):
        validate_constraint(fns.g,'g')
    if(pclass == ProblemClass.InequalityConstrained or pclass == ProblemClass.Constrained):
        validate_constraint(fns.h,'h')

# Determine the type of problem from the points involved.
def type_from_pts(pts):
    try:
        if not hasattr(pts,'x'):
            raise PeoptError('Every collection of points must include x.')
        if hasattr(pts,'y') and hasattr(pts,'z'):
            return ProblemClass.Constrained
        elif hasattr(pts,'y'):
            return ProblemClass.EqualityConstrained
        elif hasattr(pts,'z'):
            return ProblemClass.InequalityConstrained
        else:
            return ProblemClass.Unconstrained
    except PeoptError as e:
        print e.value + '\n'
        raise

# Insures that we have enough points for a finite difference test
def validate_fd_pts(pts):
    # Figure out the type of optimization problem
    pclass=type_from_pts(pts)

    try:
        # Determine what points are required to run the finite difference test
        if pclass==ProblemClass.Unconstrained:
            req_pts = ['x','dx','dxx']
        elif pclass==ProblemClass.EqualityConstrained:
            req_pts = ['x','dx','dxx','y','dy']
        elif pclass==ProblemClass.InequalityConstrained:
            req_pts = ['x','dx','dxx','z','dz']
        else:
            req_pts = ['x','dx','dxx','y','dy','z','dz']
        
        # Check that each of the points is present in the collection
        for pt in req_pts:
            if not hasattr(pts,pt):
                raise PeoptError('Missing the ' + pt + ' point.')

    except PeoptError as e:
        msg = reduce(lambda x,y:x + '    ' + y + ' = ...\n',
            req_pts,'class pts :\n')
        print('Error in the collection of points for the finite difference ' +
            'test.\n\n' + e.value + '\n\n' +
            'We require the collection of points to be a class of the form:\n'+
            msg)
        raise

# Run a series of finite difference checks and other diagnostics on the
# functions
def diagnostics(vs,fns,pts):
    # First, validate the vector space, functions, and points
    validate_vector_spaces(vs)
    validate_functions(fns)
    validate_fd_pts(pts)

    # Next, get the problem class determined by each of these structions
    vs_class = type_from_vs(vs)
    fns_class = type_from_fns(fns)
    pts_class = type_from_pts(pts)

    # Throw an exception if the problem class is not all the same
    if not(vs_class == fns_class and fns_class == pts_class):
        msg = 'The problem class differs between the vector space, ' + \
            'functions, and points.'
        print msg + '\n'
        raise PeoptError(msg)

    # If they are all the same, run the diagnostics 
    else:
        ret = libpeopt.pypeopt(vs_class,vs,fns,pts,None)
        
        # Check for an exception
        if ret[1]!=None:
            raise PeoptError(ret[1]);

# Solve the optimization problem
def getMin(vs,fns,pts,fname):
    # First, validate the vector space, functions, and points
    validate_vector_spaces(vs)
    validate_functions(fns)

    # Next, get the problem class determined by each of these structions
    vs_class = type_from_vs(vs)
    fns_class = type_from_fns(fns)
    pts_class = type_from_pts(pts)

    # Throw an exception if the problem class is not all the same
    if not(vs_class == fns_class and fns_class == pts_class):
        msg = 'The problem class differs between the vector space, ' + \
            'functions, and points.'
        print msg + '\n'
        raise PeoptError(msg)

    # Throw an exception if fname is not a string
    if not isinstance(fname,str):
        msg = 'Parameter number 4 to the solve function must be a string'
        print msg + '\n'
        raise PeoptError(msg)

    # If they are all the same, run the solver
    else:
        ret = libpeopt.pypeopt(vs_class,vs,fns,pts,fname)
        
        # Check for an exception
        if ret[1]!=None:
            raise PeoptError(ret[1]);
        else:
            return ret[0]

# Setup some documentation
diagnostics_summary="""
SUMMARY
Function diagnostics for peopt.

This function performs a number of diagnostics on the functions that
peopt uses for optimization.  On the objective function f, it performs a

Finite difference test on grad(f) using f,
Finite difference test on hess(f) using grad(f),
Symmetry test on hess(f) that verifies <hess(f)dx,dxx>=<dx,hess(f)dxx>.

On the constraints g and h, it performs a

Finite difference test on g'(x) using g,
Symmetry test that verifies <g'(x)dx,dy>=<dx,g'(x)*dy>,
Finite difference test on (g''(x)dx)*dy using g'(x)*.
"""

getMin_summary="""
SUMMARY
Minimization using the optimization algorithms inside of peopt.

This function solve problems of the form

Unconstrained:          min f(x),
Equality constrained:   min f(x) st g(x) = 0,
Inequality constrained: min f(x) st h(x) >= 0,
Constrained:            min f(x) st g(x) = 0, h(x) >= 0. 

where >= denotes general symmetric cone constraint.  This includes
linear inequalities.  In order to minimize these problems, we use a
variety of gradient based algorithms, which includes second-order
methods based on Newton's method, but also includes first-order methods
such as BFGS, and SR1 as well as simple methods such as steepest
descent.  For the inequality constraints, we use a primal-dual interior
point method.
"""

both_args="""
ARGUMENTS
vs: A class instance of the form
    class vs:
        X =  # Must be equality to an instance of a vector space. 
        Y =  # Required only if the problem is equality constrained.
        Z =  # Required only if the problem is inequality constrained.

    class SomeVectorSpace:       # Includes X, Y, and Z.
        def copy(self,x):        # copy <- x (Shallow.  No memory 
                                 # allocation.)
        def scal(self,alpha,x):  # scal <- alpha * x 
        def zero(self,x):        # zero <- 0 
        def axpy(self,alpha,x,y):# axpy <- alpha * x + y
        def innr(self,x,y):      # innr <- <x,y>

        # The rest are required for vector spaces used in Z 
        def prod(self,x,y):  # Jordan product, prod <- x o y
        def id(self,x):      # Identity element, id <- e such that
                             # x o e = x
        def linv(self,x,y):  # Jordan product inverse, linv<-inv(L(x)) y
                             # where L(x) y = x o y
        def barr(self,x):    # Barrier function, barr <- barr(x) where
                             # x o grad barr(x) = e 
        def srch(self,x,y):  # Line search, srch<-argmax {alpha in Real
                             # >= 0 : alpha x + y >= 0} where y > 0.
                             # If the argmax is infinity, then return
                             # -1.0.

fns: A class instance of the form
    class fns:    
        f = Objective()
        g = EqualityConstraint()   # Optional
        h = InequalityConstraint() # Optional

    class Objective:
        def eval(self,x):       # eval <- f(x)
        def grad(self,x):       # grad <- grad(f)(x)
        def hessvec(self,x,dx): # hessvec <- hess(f)(x)dx

    # Applies to both EqualityConstraint and InequalityConstraint.  For
    # InequalityConstraint, the function must be affine.  This means
    # that (h''(x)dx)*dy must be 0.
    class Constraint:
        def eval(self,x):       # eval <- g(x) 
        def p(self,x,dx):       # p <- g'(x)dx
        def ps(self,x,dy):      # ps <- g'(x)*dy
        def pps(self,x,dx,dy):  # pps <- (g''(x)dx)*dy
"""

diagnostics_pts = """ 
pts: A class instance of the form
    class pts:
        x= 
        dx= 
        dxx= 
        y=    # Required only if the problem is equality constrained. 
        dy=   # Required only if the problem is equality constrained. 
        z=    # Required only if the problem is inequality constrained. 
        dz=   # Required only if the problem is inequality constrained.
"""

getMin_pts = """ 
pts: A class instance of the form
    class pts:
        x= 
        y=    # Required only if the problem is equality constrained. 
        z=    # Required only if the problem is inequality constrained. 
"""

getMin_fname= """
fname: A string denoting the location and name of the .peopt parameter
       file.
"""

both_args_note = """
Note: This function automatically determines whether the problem is
unconstrained, equality constrained, inequality constrained, or fully
constrained from its arguments.  Specifically, if the optional class
members above are missing, then the function assumes that the problem
class is restricted.  For example, if we define the class fns above as
    class fns:    
        f = Objective()
        g = EqualityConstraint()
then we assume that the problem is equality constrained.  Alternatively,
if we define the class vs above as
    class vs:
        X = MyVectorSpace()
        Z = MyVectorSpace()
then we assume that the problem is inequality constrained.  The classes
vs, fns, and pts must have a compatible structure otherwise the
exception PeoptError will be thrown.  In order to be compatible, see the
comments in the definition of the classes above. 
"""

diagnostics_return="""
RETURN VALUE
None
"""

getMin_return="""
RETURN VALUE
sol: A triple of the form (x,y,z) that contains the solution of the
     optimization as well as the optimal Lagrange multipliers (dual
     variables).  In this triple, x denotes the optimal solution, y
     denotes the optimal dual variable in an equality constrained
     problem, and z denotes the optimal dual variable in an inequality
     constrained problem.  When the problem is not equality constrained,
     y is set None.  Similarly, when the problem is not inequality
     constrained, z is set to None.  Finally, the returned variables x,
     y, and z are guaranteed to have the same struture as the points 
     pts.x, pts.y, and pts.z.
"""
both_exceptions="""
EXCEPTIONS
As mentioned above, if the classes vs, pts, and fns are not compatible,
the exception PeoptError will be thrown.  In addition, if one of the 
functions defined in the vector space, the objective, or the constraints
fails, the the exception PeoptError will also be thrown with a message
indicating where.  Note, sometimes this message is deceptive.  For 
example, if we forget to return a value in the copy function inside the
vector space, then we may see an error thrown by the gradient.  This is
problematic, but difficult to detect.  Due to the differences in the
type systems between C++ and Python, there's not a good way to tell that
a value returned by, say copy, is invalid until we use it.  As such, 
systematic testing of all of these functions is recommended.  In
addition, we do not report any exceptions thrown by these functions.
If an exception is thrown inside of the vector space or one of the 
function definitions, it is treated as though these functions returned 
an invalid value and a PeoptError exception will be thrown.
"""

diagnostics.__doc__ = diagnostics_summary + both_args + diagnostics_pts + \
    both_args_note + diagnostics_return + both_exceptions
getMin.__doc__ = getMin_summary + both_args + getMin_pts + getMin_fname + \
    both_args_note + getMin_return + both_exceptions
