# Set the exposed functions
__all__ = ['peopt_fd']

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
