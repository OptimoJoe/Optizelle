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
    history_reset = Optizelle.createNatProperty(
        "history_reset", 
            "Number of failed iterations before we reset the "
            "history for quasi-Newton methods")
    iter = Optizelle.createNatProperty(
        "iter",
        "Current iteration")
    iter_max = Optizelle.createNatProperty(
        "iter_max",
        "Maximum number of optimization iterations")
    opt_stop = Optizelle.createEnumProperty(
        "opt_stop",
        Optizelle.StoppingCondition,
        "Why we've stopped the optimization")
    krylov_iter = Optizelle.createNatProperty(
        "krylov_iter",
        "Current number of Krylov iterations taken")
    krylov_iter_max = Optizelle.createNatProperty(
        "krylov_iter_max",
        "Maximum number of iterations in the Krylov method")
    krylov_iter_total = Optizelle.createNatProperty(
        "krylov_iter_total",
        "Total number of Krylov iterations taken")
    krylov_orthog_max = Optizelle.createNatProperty(
        "krylov_orthog_max",
        ("The maximum number of vectors we orthogonalize "
        "against in the Krylov method.  For something like "
        "CG, this is 1."))
    krylov_stop = Optizelle.createEnumProperty(
        "krylov_stop",
        Optizelle.KrylovStop,
        "Why the Krylov method was last stopped")
    krylov_rel_err = Optizelle.createFloatProperty(
        "krylov_rel_err",
        "Relative error in the Krylov method")
    eps_krylov = Optizelle.createFloatProperty(
        "eps_krylov",
        "Stopping tolerance for the Krylov method")
    krylov_solver = Optizelle.createEnumProperty(
        "krylov_solver",
        Optizelle.KrylovSolverTruncated, 
        "Truncated Krylov solver")
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
    failed_safeguard_max = Optizelle.createNatProperty(
        "failed_safeguard_max",
        "Number of failed safe-guard steps before quitting the method")
    failed_safeguard = Optizelle.createNatProperty(
        "failed_safeguard",
        "Number of failed safeguard steps during the last iteration")
    failed_safeguard_total = Optizelle.createNatProperty(
        "failed_safeguard_total",
        "Total number of failed safeguard steps")
    alpha_x = Optizelle.createFloatProperty(
        "alpha_x",
        ("Amount we truncate dx in order to maintain feasibility "
        " with respect to the safeguard, which probably relates to "
        "the inequality constraint"))
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
    rejected_trustregion = Optizelle.createNatProperty(
        "rejected_trustregion",
        "Number of rejected trust-region steps")
    alpha0 = Optizelle.createFloatProperty(
        "alpha0",
        "Base line-search step length")
    alpha = Optizelle.createFloatProperty(
        "alpha",
        "Actual line-search step length")
    c1 = Optizelle.createFloatProperty(
        "c1",
        "Parameter that helps govern the sufficient decrease")
    linesearch_iter = Optizelle.createNatProperty(
        "linesearch_iter",
        "Current number of iterations used in the line-search")
    linesearch_iter_max = Optizelle.createNatProperty(
        "linesearch_iter_max",
        "Maximum number of iterations used in the line-search")
    linesearch_iter_total = Optizelle.createNatProperty(
        "linesearch_iter_total",
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

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type Unconstrained.State.t."
            % (name))
