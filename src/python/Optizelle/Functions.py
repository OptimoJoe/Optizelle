# Different kinds of functions that Optizelle uses

#---ScalarValuedFunction0---
class ScalarValuedFunction(object):
    """A simple scalar valued function interface, f : X -> R"""

    def _err(self,fn):
        """Produces an error message for an undefined function"""
        raise Exception.t("%s function is not defined in a " % (fn) +
            "ScalarValuedFunction")

    def eval(self,x):
        """<- f(x)"""
        _err(self,"eval")

    def grad(self,x,grad):
        """<- grad f(x)"""
        _err(self,"grad")

    def hessvec(self,x,dx,H_dx):
        """<- hess f(x) dx"""
        _err(self,"grad")
#---ScalarValuedFunction1---

#---VectorValuedFunction0---
class VectorValuedFunction(object):
    """A vector valued function interface, f : X -> Y"""

    def _err(self,fn):
        """Produces an error message for an undefined function"""
        raise Exception.t("%s function is not defined in a " % (fn) +
            "VectorValuedFunction")

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
#---VectorValuedFunction1---

#---Operator0---
class Operator(object):
    """A linear operator specification, A : X->Y"""

    def _err(self,fn):
        """Produces an error message for an undefined function"""
        raise Exception.t("%s function is not defined in an " % (fn) +
            "Operator")

    def eval(self,state,x,y):
        """y <- A(x)"""
        _err(self,"eval")
#---Operator1---

#---StateManipulator0---
class StateManipulator(object):
    """A function that has free reign to manipulate or analyze the state"""
    def eval(self,fns,state,loc):
        """Application"""
        pass
#---StateManipulator1---
