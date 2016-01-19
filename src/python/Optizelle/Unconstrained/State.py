__all__ = [
    "t"
]

import Optizelle
import Optizelle.Utility
import copy

def allocateVectors(self,X,x):
    """Allocates memory for the state vectors"""
    self.x=X.init(x)
    self.grad=X.init(x)
    self.dx=X.init(x)
    self.x_old=X.init(x)
    self.grad_old=X.init(x)
    self.dx_old=X.init(x)
        

class t(object):
    """Internal state of the optimization"""

    def __init__(self,X,msg,x):
        """Constructor"""

        # Check our arguments
        Optizelle.checkVectorSpace("X",X)
        Optizelle.checkMessaging("msg",msg)

        # Allocate memory for our vectors
        allocateVectors(self,X,x)

        # Create the state
        Optizelle.Utility.UnconstrainedStateCreate(self,X,msg,x)

    # Create all of the properties
    eps_grad = Optizelle.createFloatProperty(
        "eps_grad",
        "Tolerance for the gradient stopping condition")
    eps_dx = Optizelle.createFloatProperty(
        "eps_dx",
        "Tolerance for the step length stopping criteria")
    algorithm_class = Optizelle.createEnumProperty(
        "algorihm_class",
        Optizelle.AlgorithmClass,
        "Algorithm class")
    stored_history = Optizelle.createNatProperty(
        "stored_history", 
        "Number of control objects to store in a quasi-Newton method")
    iter = Optizelle.createNatProperty(
        "iter",
        "Current iteration")
    iter_max = Optizelle.createNatProperty(
        "iter_max",
        "Maximum number of optimization iterations")
    glob_iter = Optizelle.createNatProperty(
        "glob_iter",
        "Globalization iteration")
    glob_iter_max = Optizelle.createNatProperty(
        "glob_iter_max",
        "Maximum number of globalization iterations before we quit")
    glob_iter_total = Optizelle.createNatProperty(
        "glob_iter_total",
        "Total number of globalization iterations taken")
    opt_stop = Optizelle.createEnumProperty(
        "opt_stop",
        Optizelle.OptimizationStop,
        "Why we've stopped the optimization")
    trunc_iter = Optizelle.createNatProperty(
        "trunc_iter",
        "Current number of truncated-CG iterations taken")
    trunc_iter_max = Optizelle.createNatProperty(
        "trunc_iter_max",
        "Maximum number of iterations used by truncated CG")
    trunc_iter_total = Optizelle.createNatProperty(
        "trunc_iter_total",
        "Total number of truncated-CG iterations taken")
    trunc_orthog_storage_max = Optizelle.createNatProperty(
        "trunc_orthog_storage_max",
        "Maximum number of vectors we orthogonalize against in truncated CG")
    trunc_orthog_iter_max = Optizelle.createNatProperty(
        "trunc_orthog_iter_max",
        "Maximum number of orthogonalization iterations in truncated CG")
    trunc_stop = Optizelle.createEnumProperty(
        "trunc_stop",
        Optizelle.TruncatedStop,
        "Why truncated CG was last stopped")
    trunc_err = Optizelle.createFloatProperty(
        "trunc_err",
        "Relative error in truncated CG")
    eps_trunc = Optizelle.createFloatProperty(
        "eps_trunc",
        "Stopping tolerance for truncated CG")
    algorithm_class = Optizelle.createEnumProperty(
        "algorithm_class",
        Optizelle.AlgorithmClass,
        "Algorithm class") 
    PH_type = Optizelle.createEnumProperty(
        "PH_type",
        Optizelle.Operators,
        "Preconditioner for the Hessian") 
    H_type = Optizelle.createEnumProperty(
        "H_type",
        Optizelle.Operators,
        "Hessian approximation") 
    norm_gradtyp = Optizelle.createFloatProperty(
        "norm_gradtyp",
        "Norm of a typical tradient") 
    norm_dxtyp = Optizelle.createFloatProperty(
        "norm_dxtyp",
        "Norm of a typical trial step") 
    x = Optizelle.createVectorProperty(
        "x",
        "Optimization variable") 
    grad = Optizelle.createVectorProperty(
        "grad",
        ("Gradient, possibly of the objective, possibly of the  Lagrangian.  "
        "It depends on the context."))
    dx = Optizelle.createVectorProperty(
        "dx",
        "Trial step") 
    x_old = Optizelle.createVectorProperty(
        "x_old",
        "Old optimization variable") 
    grad_old = Optizelle.createVectorProperty(
        "grad_old",
        "Old gradient") 
    dx_old = Optizelle.createVectorProperty(
        "dx_old",
        "Old trial step") 
    oldY = Optizelle.createVectorListProperty(
        "oldY",
        "Difference in prior gradients")
    oldS = Optizelle.createVectorListProperty(
        "oldS",
        "Difference in prior steps")
    f_x = Optizelle.createFloatProperty(
        "f_x",
        "Current value of the objective function") 
    f_xpdx = Optizelle.createFloatProperty(
        "f_xpdx",
        "Objective function at the trial step") 
    msg_level = Optizelle.createNatProperty(
        "msg_level",
        "Messaging level")
    safeguard_failed_max = Optizelle.createNatProperty(
        "safeguard_failed_max",
        "Number of failed safe-guard steps before quitting the method")
    safeguard_failed = Optizelle.createNatProperty(
        "safeguard_failed",
        "Number of failed safeguard steps during the last iteration")
    safeguard_failed_total = Optizelle.createNatProperty(
        "safeguard_failed_total",
        "Total number of failed safeguard steps")
    alpha_x = Optizelle.createFloatProperty(
        "alpha_x",
        ("Amount we truncate dx in order to maintain feasibility "
        " with respect to the safeguard, which probably relates to "
        "the inequality constraint"))
    alpha_x_qn = Optizelle.createFloatProperty(
        "alpha_x_qn",
        ("Amount we truncate dx_n in order to maintain feasibility "
        "with respect to the safeguard, which probably relates to "
        "the inequailty constraint"))
    delta = Optizelle.createFloatProperty(
        "delta",
        "Trust region radius")
    eta1 = Optizelle.createFloatProperty(
        "eta1",
        "Trust-region parameter for checking whether a step has been accepted")
    eta2 = Optizelle.createFloatProperty(
        "eta2",
        ("Trust-region parameter for checking whether we enlarge the "
        "trust-region radius"))
    ared = Optizelle.createFloatProperty(
        "ared",
        "Actual reduction")
    pred = Optizelle.createFloatProperty(
        "pred",
        "Predicted reduction")
    alpha0 = Optizelle.createFloatProperty(
        "alpha0",
        "Base line-search step length")
    alpha = Optizelle.createFloatProperty(
        "alpha",
        "Actual line-search step length")
    c1 = Optizelle.createFloatProperty(
        "c1",
        "Parameter that helps govern the sufficient decrease")
    ls_iter = Optizelle.createNatProperty(
        "ls_iter",
        "Current number of iterations used in the line-search")
    ls_iter_max = Optizelle.createNatProperty(
        "ls_iter_max",
        "Maximum number of iterations used in the line-search")
    ls_iter_total = Optizelle.createNatProperty(
        "ls_iter_total",
        "Total number of line-search iterations computed")
    eps_ls = Optizelle.createFloatProperty(
        "eps_ls",
        "Stopping tolerance for the line-search")
    dir = Optizelle.createEnumProperty(
        "dir",
        Optizelle.LineSearchDirection,
        "Search direction type")
    kind = Optizelle.createEnumProperty(
        "kind",
        Optizelle.LineSearchKind,
        "Type of line-search")
    f_diag = Optizelle.createEnumProperty(
        "f_diag",
        Optizelle.FunctionDiagnostics,
        "Function diagnostics on f")
    L_diag = Optizelle.createEnumProperty(
        "L_diag",
        Optizelle.FunctionDiagnostics,
        "Function diagnostics on the Lagrangian")
    x_diag = Optizelle.createEnumProperty(
        "x_diag",
        Optizelle.VectorSpaceDiagnostics,
        "Vector space diagnostics on X")
    dscheme = Optizelle.createEnumProperty(
        "dscheme",
        Optizelle.DiagnosticScheme,
        "Diagnostic scheme")
    eps_kind = Optizelle.createEnumProperty(
        "eps_kind",
        Optizelle.ToleranceKind,
        "Kind of stopping tolerance")

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Unconstrained.State.t."
            % (name))
