# This helps with the messaging function
from __future__ import print_function
import sys
import inspect
import numpy
import math 
import copy
import random

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
    "FunctionDiagnostics",
    "DiagnosticScheme",

    "ScalarValuedFunction",
    "VectorValuedFunction",
    "Operator",
    "StateManipulator",

    "Exception",

    "Messaging",

    "Rm"
]
__doc__ = "Optizelle optimization library."

class Exception(Exception):
    """General purpose exception for all of Optizelle's errors"""
    pass

class EnumeratedType(object):
    """A generic enumerated type"""
    @classmethod
    def to_string(cls,i):
        """Converts the enumerated type into a string"""
        return filter(lambda (name,value):value==i,
            cls.__dict__.items())[0][0]
        
class KrylovStop(EnumeratedType):
    """Reasons we stop the Krylov method"""
    NegativeCurvature, \
    RelativeErrorSmall, \
    MaxItersExceeded, \
    TrustRegionViolated, \
    Instability, \
    InvalidTrustRegionCenter \
     = range(6)

class AlgorithmClass(EnumeratedType):
    """Which algorithm class do we use"""
    TrustRegion, \
    LineSearch, \
    UserDefined \
     = range(3)

class StoppingCondition(EnumeratedType):
    """Reasons why we stop the algorithm"""
    NotConverged, \
    RelativeGradientSmall, \
    RelativeStepSmall, \
    MaxItersExceeded, \
    InteriorPointInstability, \
    UserDefined \
    = range(6)

class Operators(EnumeratedType):
    """Various operators for both Hessian approximations and preconditioners"""
    Identity, \
    ScaledIdentity, \
    BFGS, \
    InvBFGS, \
    SR1, \
    InvSR1, \
    UserDefined \
    = range(7)
    
class LineSearchDirection(EnumeratedType):
    """Different kinds of search directions"""
    SteepestDescent, \
    FletcherReeves, \
    PolakRibiere, \
    HestenesStiefel, \
    BFGS, \
    NewtonCG \
    = range(6)
   
class LineSearchKind(EnumeratedType):
    """Different sorts of line searches"""
    Brents, \
    GoldenSection, \
    BackTracking, \
    TwoPointA, \
    TwoPointB \
    = range(5)
    
class OptimizationLocation(EnumeratedType):
    """Different points in the optimization algorithm"""
    BeginningOfOptimization, \
    BeforeInitialFuncAndGrad, \
    AfterInitialFuncAndGrad, \
    BeforeOptimizationLoop, \
    BeginningOfOptimizationLoop, \
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
    = range(20)

class ProblemClass(EnumeratedType):
    """Different problem classes"""
    Unconstrained, \
    EqualityConstrained, \
    InequalityConstrained, \
    Constrained \
    = range(4)
    
class KrylovSolverTruncated(EnumeratedType):
    """Different truncated Krylov solvers"""
    ConjugateDirection, \
    MINRES \
    = range(2)

class InteriorPointMethod(EnumeratedType):
    """Different kinds of interior point methods"""
    PrimalDual, \
    PrimalDualLinked, \
    LogBarrier \
    = range(3)
    
class CentralityStrategy(EnumeratedType):
    """Different schemes for adjusting the interior point centrality"""
    Constant, \
    StairStep, \
    PredictorCorrector \
    = range(3)

class FunctionDiagnostics(EnumeratedType):
    """Different function diagnostics on the optimization functions""" 
    NoDiagnostics, \
    FirstOrder, \
    SecondOrder \
    = range(3)

class DiagnosticScheme(EnumeratedType):
    """When and how often we compute our intrusive diagnostics"""
    Never, \
    DiagnosticsOnly, \
    EveryIteration \
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
        raise TypeError("The %s member must be an enumerated type (natural.)"
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

def checkStaticMethod(vsname,name,value):
    """Check that we have a method"""
    if not (
        hasattr(value,name) and isinstance(value.__dict__[name],staticmethod)):
        raise TypeError("The %s member is required as a static member in " ^
            "the vector space %s." % (name,vsname))

def checkVectorSpace(vsname,value):
    """Check that we have a valid-vector space"""

    # Define all the functions we care about
    fns=["init","copy","scal","zero","axpy","innr","rand"]

    # Now, check each of these
    map(lambda name:checkStaticMethod(vsname,name,value),fns) 

def checkEuclidean(vsname,value):
    """Check that we have a valid Euclidean-Jordan algebra"""

    # Check that we have a valid vector space
    checkVectorSpace(vsname,value)
    
    # Define all the new functions we care about
    fns=["prod","id","linv","barr","srch","symm"]

    # Now, check each of these
    map(lambda name:checkStaticMethod(vsname,name,value),fns) 

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

def checkVectors(name,value):
    """Check that we have a list of restart vectors"""
    if not issubclass(type(value),list):
        raise TypeError("The %s argument must be a list." % (name))
    map(lambda (i,x):checkString("%s[%d][0]" % (name,i),x[0]),
        enumerate(value))

def checkReals(name,value):
    """Check that we have a list of restart reals"""
    if not issubclass(type(value),list):
        raise TypeError("The %s argument must be a list." % (name))
    map(lambda (i,x):checkString("%s[%d][0]" % (name,i),x[0]),
        enumerate(value))
    map(lambda (i,x):checkFloat("%s[%d][1]" % (name,i),x[1]),
        enumerate(value))

def checkNaturals(name,value):
    """Check that we have a list of restart naturals"""
    if not issubclass(type(value),list):
        raise TypeError("The %s argument must be a list." % (name))
    map(lambda (i,x):checkString("%s[%d][0]" % (name,i),x[0]),
        enumerate(value))
    map(lambda (i,x):checkNatural("%s[%d][1]" % (name,i),x[1]),
        enumerate(value))

def checkParams(name,value):
    """Check that we have a list of restart parameters"""
    if not issubclass(type(value),list):
        raise TypeError("The %s argument must be a list." % (name))
    map(lambda (i,x):checkString("%s[%d][0]" % (name,i),x[0]),
        enumerate(value))
    map(lambda (i,x):checkString("%s[%d][1]" % (name,i),x[1]),
        enumerate(value))

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
    """Vector space for the nonnegative orthant.  For basic vectors in R^m, use this."""

    @staticmethod
    def init(x):
        """Memory allocation and size setting"""
        return copy.deepcopy(x) 

    @staticmethod
    def copy(x,y):
        """y <- x (Shallow.  No memory allocation.)"""
        numpy.copyto(y,x) 

    @staticmethod
    def scal(alpha,x):
        """x <- alpha * x"""
        x.__imul__(alpha)

    @staticmethod
    def zero(x):
        """x <- 0"""
        x.fill(0.)

    @staticmethod
    def axpy(alpha,x,y):
        """y <- alpha * x + y"""
        y.__iadd__(alpha*x)

    @staticmethod
    def innr(x,y):
        """<- <x,y>"""
        return numpy.inner(x,y) 

    @staticmethod
    def rand(x):
        """x <- random"""
        numpy.copyto(x,map(lambda x:random.normalvariate(0.,1.),x))

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

class RestartPackage(tuple):
    """Holds restart information"""
    def __new__ (cls):
        return super(RestartPackage,cls).__new__(cls,tuple([[],[]]))
    def __init__(self):
        super(RestartPackage,self).__init__(tuple([[],[]]))
