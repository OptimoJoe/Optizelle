# Modules required for the functions here
import sys
import inspect
import numpy
import math
import copy
import random
import collections
import functools

# Import the Optizelle pieces, which actually depend on this module
from Optizelle.Enumerated import *
from Optizelle.Functions import *

import Optizelle.Unconstrained.State
import Optizelle.Unconstrained.Functions
import Optizelle.Unconstrained.Algorithms
import Optizelle.Unconstrained.Restart

import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions
import Optizelle.EqualityConstrained.Algorithms
import Optizelle.EqualityConstrained.Restart

import Optizelle.InequalityConstrained.State
import Optizelle.InequalityConstrained.Functions
import Optizelle.InequalityConstrained.Algorithms
import Optizelle.InequalityConstrained.Restart

import Optizelle.Constrained.State
import Optizelle.Constrained.Functions
import Optizelle.Constrained.Algorithms
import Optizelle.Constrained.Restart

import Optizelle.json.Serialization
import Optizelle.json.Unconstrained
import Optizelle.json.EqualityConstrained
import Optizelle.json.InequalityConstrained
import Optizelle.json.Constrained

import Optizelle.Messaging
import Optizelle.Exception

__all__ = [
    "Unconstrained",
    "InequalityConstrained",
    "EqualityConstrained",
    "Constrained",
    "Utility"

    "TruncatedStop",
    "AlgorithmClass",
    "OptimizationStop",
    "Operators",
    "LineSearchDirection",
    "LineSearchKind",
    "OptimizationLocation",
    "ProblemClass",
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
__doc__ = "Optizelle optimization library"

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
        numpy.copyto(x,numpy.vectorize(lambda x:random.normalvariate(0.,1.))(x))

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
        if (x>0).all():
            return functools.reduce(lambda x,y:x+math.log(y),x,0.)
        else:
            return float("nan")

    @staticmethod
    def srch(x,y):
        """Line search, <- argmax {alpha \in Real >= 0 : alpha x + y >= 0} where y > 0"""
        alpha = float("inf")
        for i in range(0,len(x)):
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
