# This helps with the messaging function
from __future__ import print_function
import sys
import inspect
import numpy
import math 
import copy

__all__ = [
    "Unconstrained",
    "InequalityConstrained",
    "EqualityConstrained",
    "Constrained",
    "Utility"

    "KrylovStop",
    "AlgorithmClass",
    "StoppingCondition",
    "Operators",
    "LineSearchDirection",
    "LineSearchKind",
    "OptimizationLocation",
    "ProblemClass",
    "KrylovSolverTruncated",
    "InteriorPointMethod",
    "CentralityStrategy",

    "ScalarValuedFunction",
    "VectorValuedFunction",
    "Operator",

    "Exception",

    "Messaging",

    "Rm"
]
__doc__ = "Optizelle optimization library."

class Exception(Exception):
    """General purpose exception for all of Optizelle's errors"""
    pass

class KrylovStop(object):
    """Reasons we stop the Krylov method"""
    NegativeCurvature, \
    RelativeErrorSmall, \
    MaxItersExceeded, \
    TrustRegionViolated, \
    Instability, \
    InvalidTrustRegionCenter \
     = range(6)

class AlgorithmClass(object):
    """Which algorithm class do we use"""
    TrustRegion, \
    LineSearch, \
    UserDefined \
     = range(3)

class StoppingCondition(object):
    """Reasons why we stop the algorithm"""
    NotConverged, \
    RelativeGradientSmall, \
    RelativeStepSmall, \
    MaxItersExceeded, \
    InteriorPointInstability, \
    UserDefined \
    = range(6)

class Operators(object):
    """Various operators for both Hessian approximations and preconditioners"""
    Identity, \
    ScaledIdentity, \
    BFGS, \
    InvBFGS, \
    SR1, \
    InvSR1, \
    UserDefined \
    = range(7)
    
class LineSearchDirection(object):
    """Different kinds of search directions"""
    SteepestDescent, \
    FletcherReeves, \
    PolakRibiere, \
    HestenesStiefel, \
    BFGS, \
    NewtonCG \
    = range(6)
   
class LineSearchKind(object):
    """Different sorts of line searches"""
    Brents, \
    GoldenSection, \
    BackTracking, \
    TwoPointA, \
    TwoPointB \
    = range(5)
    
class OptimizationLocation(object):
    """Different points in the optimization algorithm"""
    BeginningOfOptimization, \
    BeforeInitialFuncAndGrad, \
    AfterInitialFuncAndGrad, \
    BeforeOptimizationLoop, \
    BeforeSaveOld, \
    BeforeStep, \
    BeforeGetStep, \
    GetStep, \
    AfterStepBeforeGradient, \
    AfterGradient, \
    BeforeQuasi, \
    AfterQuasi, \
    EndOfOptimizationIteration, \
    BeforeLineSearch, \
    AfterRejectedTrustRegion, \
    AfterRejectedLineSearch, \
    BeforeActualVersusPredicted, \
    EndOfKrylovIteration, \
    EndOfOptimization \
    = range(19)

class ProblemClass(object):
    """Different problem classes"""
    Unconstrained, \
    EqualityConstrained, \
    InequalityConstrained, \
    Constrained \
    = range(4)
    
class KrylovSolverTruncated(object):
    """Different truncated Krylov solvers"""
    ConjugateDirection, \
    MINRES \
    = range(2)

class InteriorPointMethod(object):
    """Different kinds of interior point methods"""
    PrimalDual, \
    PrimalDualLinked, \
    LogBarrier \
    = range(3)
    
class CentralityStrategy(object):
    """Different schemes for adjusting the interior point centrality"""
    Constant, \
    StairStep, \
    PredictorCorrector \
    = range(3)

def checkFloat(name,value):
    """Checks that an input is a floating-point number"""
    if type(value)!=float:
        raise TypeError("The %s member must be a floating point." % name)

def checkNatural(name,value):
    """Checks that an input is a natural number"""
    if type(value)!=int or value < 0:
        raise TypeError("The %s member must be a natural number." % name)

def checkEnum(name,value):
    """Checks that an input is an enumerated type """
    if type(value)!=int or value < 0:
        raise TypeError("The %s member must be an enumerated type (natural)."
            % name)

def checkEnumRange(name,enum,value):
    """Checks that an input is in a valid enumerated range""" 
    if not value in enum.__dict__.itervalues(): 
        raise TypeError("The %s member is outside the valid enumated range."
            % name)

def checkVectorList(name,value):
    """Checks that an input is a list"""
    if not type(value)==list:
        raise TypeError("The %s member must be a list of vectors." % name)

def checkFunction(name,value):
    """Checks that an input is a function"""
    if not inspect.isfunction(value):
        raise TypeError("The %s member must be a function." % name)

def checkDelete(name):
    """Check that we don't delete something"""
    raise TypeError("Cannot delete the %s member." % name)

def checkScalarValuedFunction(name,value):
    """Check that we have a scalar-valued function"""
    if not issubclass(type(value),ScalarValuedFunction):
        raise TypeError("The %s member must be a ScalarValuedFunction." % name)

def checkVectorValuedFunction(name,value):
    """Check that we have a vector-valued function"""
    if not issubclass(type(value),VectorValuedFunction):
        raise TypeError("The %s member must be a VectorValuedFunction." % name)

def checkOperator(name,value):
    """Check that we have a linear operator"""
    if not issubclass(type(value),Operator):
        raise TypeError("The %s member must be an Operator." % name)

def checkMethod(vsname,name,value):
    """Check that we have a method"""
    if not (hasattr(value,name) and inspect.isfunction(value.__dict__[name])):
        raise TypeError("The %s member is required in the vector space %s."
            % (name,vsname))

def checkVectorSpace(vsname,value):
    """Check that we have a valid-vector space"""

    # Define all the functions we care about
    fns=["copy","scal","zero","axpy","innr"]

    # Now, check each of these
    map(lambda name:checkMethod(vsname,name,value),fns) 

def checkEuclidean(vsname,value):
    """Check that we have a valid Euclidean-Jordan algebra"""

    # Check that we have a valid vector space
    checkVectorSpace(vsname,value)
    
    # Define all the new functions we care about
    fns=["prod","id","linv","barr","srch","symm"]

    # Now, check each of these
    map(lambda name:checkMethod(vsname,name,value),fns) 

def checkMessaging(name,value):
    """Check that we have a messaging object"""
    if not issubclass(type(value),Messaging):
        raise TypeError("The %s argument must be a Messaging object." % (name))

def checkString(name,value):
    """Check that we have a string object"""
    if type(value)!=str:
        raise TypeError("The %s argument must be a string." % (name))

def checkStateManipulator(name,value):
    """Check that we have a state manipulator""" 
    if not issubclass(type(value),StateManipulator):
        raise TypeError("The %s argument must be a StateManipulator object."
            % (name))

def checkType(name,value):
    """Check that we have a type"""
    if type(value)!=type(type):
        raise TypeError("The %s argument must be a type." % (name))

def createFloatProperty(name,desc):
    """Create a floating-point property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkFloat(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createNatProperty(name,desc):
    """Create a natural number property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkNatural(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createEnumProperty(name,enum,desc):
    """Create an enumerated type property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkEnum(name,value)
        checkEnumRange(name,enum,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createFunctionProperty(name,desc):
    """Create a function property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkFunction(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createVectorProperty(name,desc):
    """Create a vector property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createVectorListProperty(name,desc):
    """Create a list of vectors property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkVectorList(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)


def createScalarValuedFunctionProperty(name,desc):
    """Create a scalar-valued function property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkScalarValuedFunction(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createVectorValuedFunctionProperty(name,desc):
    """Create a vector-valued function property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkVectorValuedFunction(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

def createOperatorProperty(name,desc):
    """Create an operator property"""
    def getter(self):
        return self.__dict__["_%s" % name] 

    def setter(self, value):
        checkOperator(name,value)
        self.__dict__["_%s" % name] = value

    def deleter(self):
        checkDelete(name)

    return property(getter,setter,deleter,desc)

class ScalarValuedFunction(object):
    """A simple scalar valued function interface, f : X -> R"""

    def _err(self,fn):
        """Produces an error message for an undefined function."""
        raise Exception("The %s function is not defined in a " % (fn) +
            "ScalarValuedFunction.")

    def eval(self,x):
        """<- f(x)"""
        _err(self,"eval")
    
    def grad(self,x,grad):
        """<- grad f(x)"""
        _err(self,"grad")
    
    def hessvec(self,x,dx,H_dx):
        """<- hess f(x) dx"""
        _err(self,"grad")

class VectorValuedFunction(object):
    """A vector valued function interface, f : X -> Y"""

    def _err(self,fn):
        """Produces an error message for an undefined function."""
        raise Exception("The %s function is not defined in a " % (fn) +
            "VectorValuedFunction.")

    def eval(self,x,y):
        """y <- f(x)"""
        _err(self,"eval")
    
    def p(self,x,dx,y):
        """y <- f'(x)dx"""
        _err(self,"p")
    
    def ps(self,x,dx,z):
        """z <- f'(x)dx"""
        _err(self,"ps")
    
    def pps(self,x,dx,dy,z):
        """z <- (f''(x)dx)*dy"""
        _err(self,"pps")
    
class Operator(object):
    """A linear operator specification, A : X->Y"""
    
    def _err(self,fn):
        """Produces an error message for an undefined function."""
        raise Exception("The %s function is not defined in an " % (fn) +
            "Operator.")

    def eval(self,state,x,y):
        """y <- A(x)"""
        _err(self,"eval")

class Messaging(object):
    """Defines how we output messages to the user"""
   
    def print(self,msg):
        """Prints out normal diagnostic information"""
        sys.stdout.write("%s\n" %(msg))

    def error(self,msg):
        """Prints out error information"""
        sys.stderr.write("%s\n" %(msg))

class StateManipulator(object):
    """A function that has free reign to manipulate or analyze the state.  This should be used cautiously."""
    def eval(self,fns,state,loc):
        """Application"""
        pass

class Rm(object):
    """ Vector space for the nonnegative orthant.  For basic vectors in R^m, use this."""

    def init(x):
        """Memory allocation and size setting"""
        return copy.deepcopy(x) 

    def copy(x,y):
        """y <- x (Shallow.  No memory allocation.)"""
        return numpy.copyto(y,x) 

    def scal(alpha,x):
        """x <- alpha * x"""
        x.__imul__(alpha)

    def zero(x):
        """x <- 0"""
        x.fill(0.)

    def axpy(alpha,x,y):
        """y <- alpha * x + y"""
        y.__iadd__(alpha*x)

    def innr(x,y):
        """<- <x,y>"""
        return numpy.inner(x,y) 

    def prod(x,y,z):
        """Jordan product, z <- x o y"""
        numpy.copyto(z,x*y)

    def id(x):
        """Identity element, x <- e such that x o e = x"""
        x.fill(1.)

    def linv(x,y,z):
        """Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y"""
        numpy.copyto(z,numpy.divide(y,x))

    def barr(x):
        """Barrier function, <- barr(x) where x o grad barr(x) = e"""
        return reduce(lambda x,y:x+math.log(y),x,0.)
        
    def srch(x,y):
        """Line search, <- argmax {alpha \in Real >= 0 : alpha x + y >= 0} where y > 0"""
        alpha = float("inf")
        for i in xrange(0,len(x)):
            if x[i] < 0:
                alpha0 = -y[i]/x[i]
                if alpha0 < alpha:
                    alpha=alpha0
        return alpha

    def symm(x):
        """Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric operator"""
        pass

class RestartPackage(tuple):
    """Holds restart information"""
    def __new__ (cls):
        return super(RestartPackage,cls).__new__(cls,tuple([[],[]]))
    def __init__(self):
        super(RestartPackage,self).__init__(tuple([[],[]]))
