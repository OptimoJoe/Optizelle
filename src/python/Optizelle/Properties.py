# Functions that create and check the elements contained in the state and
# bundle of functions

import inspect
from Optizelle.Functions import *
import pdb

__doc__ = "Optizelle support functions"

def checkFloat(name,value):
    """Checks that an input is a floating-point number"""
    if not isinstance(value,float):
        raise TypeError("%s member must be a floating point" % name)

def checkNatural(name,value):
    """Checks that an input is a natural number"""
    if not isinstance(value,int) or value < 0:
        raise TypeError("%s member must be a natural number" % name)

def checkEnum(name,value):
    """Checks that an input is an enumerated type """
    if not isinstance(value,int) or value < 0:
        raise TypeError("%s member must be an enumerated type (natural.)"
            % name)

def checkEnumRange(name,enum,value):
    """Checks that an input is in a valid enumerated range"""
    if not value in enum.__dict__.values():
        raise TypeError("%s member is outside the valid enumated range"
            % name)

def checkVectorList(name,value):
    """Checks that an input is a list"""
    if not isinstance(value,list):
        raise TypeError("%s member must be a list of vectors" % name)

def checkFunction(name,value):
    """Checks that an input is a function"""
    if not inspect.isfunction(value):
        raise TypeError("%s member must be a function" % name)

def checkDelete(name):
    """Check that we don't delete something"""
    raise TypeError("Cannot delete the %s member" % name)

def checkScalarValuedFunction(name,value):
    """Check that we have a scalar-valued function"""
    if not issubclass(type(value),ScalarValuedFunction):
        raise TypeError("%s member must be a ScalarValuedFunction" % name)

def checkVectorValuedFunction(name,value):
    """Check that we have a vector-valued function"""
    if not issubclass(type(value),VectorValuedFunction):
        raise TypeError("%s member must be a VectorValuedFunction" % name)

def checkOperator(name,value):
    """Check that we have a linear operator"""
    if not issubclass(type(value),Operator):
        raise TypeError("%s member must be an Operator" % name)

def checkStaticMethod(vsname,name,value):
    """Check that we have a method"""
    if not (
        hasattr(value,name) and isinstance(value.__dict__[name],staticmethod)):
        raise TypeError("%s member is required as a static member in " ^
            "the vector space %s" % (name,vsname))

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
    if not inspect.isfunction(value):
        raise TypeError("%s argument must be a messaging function" %(name))

def checkString(name,value):
    """Check that we have a string object"""
    if not isinstance(value,str):
        raise TypeError("%s argument must be a string" % (name))

def checkStateManipulator(name,value):
    """Check that we have a state manipulator"""
    if not issubclass(type(value),StateManipulator):
        raise TypeError("%s argument must be a StateManipulator object"
            % (name))

def checkType(name,value):
    """Check that we have a type"""
    if type(value)!=type(type):
        raise TypeError("%s argument must be a type" % (name))

def checkVectors(name,value):
    """Check that we have a list of restart vectors"""
    if not issubclass(type(value),list):
        raise TypeError("%s argument must be a list" % (name))
    map(lambda i_x:checkString("%s[%d][0]" % (name,i_x[0]),i_x[1][0]),
        enumerate(value))

def checkReals(name,value):
    """Check that we have a list of restart reals"""
    if not issubclass(type(value),list):
        raise TypeError("%s argument must be a list" % (name))
    map(lambda i_x:checkString("%s[%d][0]" % (name,i_x[0]),i_x[1][0]),
        enumerate(value))
    map(lambda i_x:checkString("%s[%d][1]" % (name,i_x[0]),i_x[1][1]),
        enumerate(value))

def checkNaturals(name,value):
    """Check that we have a list of restart naturals"""
    if not issubclass(type(value),list):
        raise TypeError("%s argument must be a list" % (name))
    map(lambda i_x:checkString("%s[%d][0]" % (name,i_x[0]),i_x[1][0]),
        enumerate(value))
    map(lambda i_x:checkString("%s[%d][1]" % (name,i_x[0]),i_x[1][1]),
        enumerate(value))

def checkParams(name,value):
    """Check that we have a list of restart parameters"""
    if not issubclass(type(value),list):
        raise TypeError("%s argument must be a list" % (name))
    map(lambda i_x:checkString("%s[%d][0]" % (name,i_x[0]),i_x[1][0]),
        enumerate(value))
    map(lambda i_x:checkString("%s[%d][1]" % (name,i_x[0]),i_x[1][1]),
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
