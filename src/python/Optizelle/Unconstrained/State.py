__all__ = [
    "t"
]

from Optizelle.Utility import *
from Optizelle.Properties import *
from Optizelle.Enumerated import *

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

    def __init__(self,X,x):
        """Constructor"""

        # Check our arguments
        checkVectorSpace("X",X)

        # Allocate memory for our vectors
        allocateVectors(self,X,x)

        # Create the state
        UnconstrainedStateCreate(self,X,x)

    # Create all of the properties
    eps_grad = createFloatProperty(
        "eps_grad",
        "Tolerance for the gradient stopping condition")
    eps_dx = createFloatProperty(
        "eps_dx",
        "Tolerance for the step length stopping criteria")
    algorithm_class = createEnumProperty(
        "algorihm_class",
        AlgorithmClass,
        "Algorithm class")
    stored_history = createNatProperty(
        "stored_history",
        "Number of control objects to store in a quasi-Newton method")
    iter = createNatProperty(
        "iter",
        "Current iteration")
    iter_max = createNatProperty(
        "iter_max",
        "Maximum number of optimization iterations")
    glob_iter = createNatProperty(
        "glob_iter",
        "Globalization iteration")
    glob_iter_max = createNatProperty(
        "glob_iter_max",
        "Maximum number of globalization iterations before we quit")
    glob_iter_total = createNatProperty(
        "glob_iter_total",
        "Total number of globalization iterations taken")
    opt_stop = createEnumProperty(
        "opt_stop",
        OptimizationStop,
        "Why we've stopped the optimization")
    trunc_iter = createNatProperty(
        "trunc_iter",
        "Current number of truncated-CG iterations taken")
    trunc_iter_max = createNatProperty(
        "trunc_iter_max",
        "Maximum number of iterations used by truncated CG")
    trunc_iter_total = createNatProperty(
        "trunc_iter_total",
        "Total number of truncated-CG iterations taken")
    trunc_orthog_storage_max = createNatProperty(
        "trunc_orthog_storage_max",
        "Maximum number of vectors we orthogonalize against in truncated CG")
    trunc_orthog_iter_max = createNatProperty(
        "trunc_orthog_iter_max",
        "Maximum number of orthogonalization iterations in truncated CG")
    trunc_stop = createEnumProperty(
        "trunc_stop",
        TruncatedStop,
        "Why truncated CG was last stopped")
    trunc_err = createFloatProperty(
        "trunc_err",
        "Relative error in truncated CG")
    eps_trunc = createFloatProperty(
        "eps_trunc",
        "Stopping tolerance for truncated CG")
    algorithm_class = createEnumProperty(
        "algorithm_class",
        AlgorithmClass,
        "Algorithm class")
    PH_type = createEnumProperty(
        "PH_type",
        Operators,
        "Preconditioner for the Hessian")
    H_type = createEnumProperty(
        "H_type",
        Operators,
        "Hessian approximation")
    norm_gradtyp = createFloatProperty(
        "norm_gradtyp",
        "Norm of a typical tradient")
    norm_dxtyp = createFloatProperty(
        "norm_dxtyp",
        "Norm of a typical trial step")
    x = createVectorProperty(
        "x",
        "Optimization variable")
    grad = createVectorProperty(
        "grad",
        ("Gradient, possibly of the objective, possibly of the  Lagrangian.  "
        "It depends on the context."))
    dx = createVectorProperty(
        "dx",
        "Trial step")
    x_old = createVectorProperty(
        "x_old",
        "Old optimization variable")
    grad_old = createVectorProperty(
        "grad_old",
        "Old gradient")
    dx_old = createVectorProperty(
        "dx_old",
        "Old trial step")
    oldY = createVectorListProperty(
        "oldY",
        "Difference in prior gradients")
    oldS = createVectorListProperty(
        "oldS",
        "Difference in prior steps")
    f_x = createFloatProperty(
        "f_x",
        "Current value of the objective function")
    f_xpdx = createFloatProperty(
        "f_xpdx",
        "Objective function at the trial step")
    msg_level = createNatProperty(
        "msg_level",
        "Messaging level")
    safeguard_failed_max = createNatProperty(
        "safeguard_failed_max",
        "Number of failed safe-guard steps before quitting the method")
    safeguard_failed = createNatProperty(
        "safeguard_failed",
        "Number of failed safeguard steps during the last iteration")
    safeguard_failed_total = createNatProperty(
        "safeguard_failed_total",
        "Total number of failed safeguard steps")
    alpha_x = createFloatProperty(
        "alpha_x",
        ("Amount we truncate dx in order to maintain feasibility "
        " with respect to the safeguard, which probably relates to "
        "the inequality constraint"))
    alpha_x_qn = createFloatProperty(
        "alpha_x_qn",
        ("Amount we truncate dx_n in order to maintain feasibility "
        "with respect to the safeguard, which probably relates to "
        "the inequailty constraint"))
    delta = createFloatProperty(
        "delta",
        "Trust region radius")
    eta1 = createFloatProperty(
        "eta1",
        "Trust-region parameter for checking whether a step has been accepted")
    eta2 = createFloatProperty(
        "eta2",
        ("Trust-region parameter for checking whether we enlarge the "
        "trust-region radius"))
    ared = createFloatProperty(
        "ared",
        "Actual reduction")
    pred = createFloatProperty(
        "pred",
        "Predicted reduction")
    alpha0 = createFloatProperty(
        "alpha0",
        "Base line-search step length")
    alpha = createFloatProperty(
        "alpha",
        "Actual line-search step length")
    c1 = createFloatProperty(
        "c1",
        "Parameter that helps govern the sufficient decrease")
    ls_iter = createNatProperty(
        "ls_iter",
        "Current number of iterations used in the line-search")
    ls_iter_max = createNatProperty(
        "ls_iter_max",
        "Maximum number of iterations used in the line-search")
    ls_iter_total = createNatProperty(
        "ls_iter_total",
        "Total number of line-search iterations computed")
    eps_ls = createFloatProperty(
        "eps_ls",
        "Stopping tolerance for the line-search")
    dir = createEnumProperty(
        "dir",
        LineSearchDirection,
        "Search direction type")
    kind = createEnumProperty(
        "kind",
        LineSearchKind,
        "Type of line-search")
    f_diag = createEnumProperty(
        "f_diag",
        FunctionDiagnostics,
        "Function diagnostics on f")
    L_diag = createEnumProperty(
        "L_diag",
        FunctionDiagnostics,
        "Function diagnostics on the Lagrangian")
    x_diag = createEnumProperty(
        "x_diag",
        VectorSpaceDiagnostics,
        "Vector space diagnostics on X")
    dscheme = createEnumProperty(
        "dscheme",
        DiagnosticScheme,
        "Diagnostic scheme")
    eps_kind = createEnumProperty(
        "eps_kind",
        ToleranceKind,
        "Kind of stopping tolerance")

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Unconstrained.State.t."
            % (name))
