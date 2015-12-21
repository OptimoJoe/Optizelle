__all__ = [
    "t"
]

import Optizelle
import Optizelle.Unconstrained.State

def allocateVectors(self,X,Y,x,y):
    """Allocates memory for the state vectors"""
    self.y=Y.init(y)
    self.dy=Y.init(y)
    self.g_x=Y.init(y)
    self.gpxdxn_p_gx=Y.init(y)
    self.gpxdxt=Y.init(y)
    self.dx_n=X.init(x)
    self.dx_ncp=X.init(x)
    self.dx_t=X.init(x)
    self.dx_t_uncorrected=X.init(x)
    self.dx_tcp_uncorrected=X.init(x)
    self.H_dxn=X.init(x)
    self.W_gradpHdxn=X.init(x)
    self.H_dxtuncorrected=X.init(x)

class t(Optizelle.Unconstrained.State.t):
    """Internal state of the optimization"""

    def __init__(self,X,Y,msg,x,y):
        """Constructor"""

        # Check our arguments
        Optizelle.checkVectorSpace("X",X)
        Optizelle.checkVectorSpace("Y",Y)
        Optizelle.checkMessaging("msg",msg)
        
        # Allocate memory for our vectors
        Optizelle.Unconstrained.State.allocateVectors(self,X,x)
        allocateVectors(self,X,Y,x,y)

        # Create the state
        Optizelle.Utility.EqualityConstrainedStateCreate(self,X,Y,msg,x,y)

    # Create all of the properties
    y = Optizelle.createVectorProperty(
        "y",
        "Equality multiplier (dual variable or Lagrange multiplier)")
    dy = Optizelle.createVectorProperty(
        "dy",
        "Step in the equality multiplier")
    zeta = Optizelle.createFloatProperty(
        "zeta",
        "The fraction of the total trust-region used for the quasi-norm step")
    eta0 = Optizelle.createFloatProperty(
        "eta0",
        ("Trust-region parameter that bounds the error in the predicted "
        "reduction"))
    rho = Optizelle.createFloatProperty(
        "rho",
        "Penalty parameter for the augmented-Lagrangian")
    rho_old = Optizelle.createFloatProperty(
        "rho_old",
        "Penalty parameter from the last iteration")
    rho_bar = Optizelle.createFloatProperty(
        "rho_bar",
        "Fixed increase in the penalty parameter")
    eps_constr = Optizelle.createFloatProperty(
        "eps_constr",
        "Stopping tolerance for the norm of the constraints")
    xi_qn = Optizelle.createFloatProperty(
        "xi_qn",
        "Inexactness tolerance for the quasi-Newton step")
    xi_pg = Optizelle.createFloatProperty(
        "xi_pg",
        "Inexactness tolerance for the projection of the gradient")
    xi_proj = Optizelle.createFloatProperty(
        "xi_proj",
        "Inexactness tolerance for the null-space projection")
    xi_tang = Optizelle.createFloatProperty(
        "xi_tang",
        "Inexactness tolerance for the tangential step")
    xi_lmh = Optizelle.createFloatProperty(
        "xi_lmh",
        "Inexactness tolerance for the equality multiplier")
    def xi_all(self,value):
        """Sets all the inexactness tolerances: xi_qn, xi_pg, xi_proj, xi_tang, and xi_lmh"""
        self.xi_qn=value
        self.xi_pg=value
        self.xi_proj=value
        self.xi_tang=value
        self.xi_lmh=value
    xi_lmg = Optizelle.createFloatProperty(
        "xi_lmg",
        "Absolute tolerance on the residual of the equality multiplier solve")
    xi_4 = Optizelle.createFloatProperty(
        "xi_4",
        ("Tolerance for how much error is acceptable after computing the "
        "tangential step given the result from the tangential subproblem"))
    rpred = Optizelle.createFloatProperty(
        "rpred",
        "Residual term in the predicted reduction")
    PSchur_left_type = Optizelle.createEnumProperty(
        "PSchur_left_type",
        Optizelle.Operators,
        "Left preconditioner for the augmented system")
    PSchur_right_type = Optizelle.createEnumProperty(
        "PSchur_right_type",
        Optizelle.Operators,
        "Right preconditioner for the augmented system")
    augsys_iter_max = Optizelle.createNatProperty(
        "augsys_iter_max",
        "Maximum number of iterations used when solving the augmented system")
    augsys_rst_freq = Optizelle.createNatProperty(
        "augsys_rst_freq",
        ("How often we restart the augmented system solve"))
    augsys_qn_iter = Optizelle.createNatProperty(
        "augsys_qn_iter",
        ("Number of augmented system solve iterations used on the quasi-normal "
        "step"))
    augsys_pg_iter = Optizelle.createNatProperty(
        "augsys_pg_iter",
        ("Number of augmented system solve iterations used when projecting the "
        "gradient prior to the tangential subproblem"))
    augsys_proj_iter = Optizelle.createNatProperty(
        "augsys_proj_iter",
        ("Number of augmented system solve iterations used in the nullspace "
        "projection inside the tangential subproblem"))
    augsys_tang_iter = Optizelle.createNatProperty(
        "augsys_tang_iter",
        ("Number of augmented system solve iterations used in the tangential "
        "step"))
    augsys_lmh_iter = Optizelle.createNatProperty(
        "augsys_lmh_iter",
        ("Number of augmented system solve iterations used in the equality "
        "multiplier solve"))
    augsys_qn_iter_total = Optizelle.createNatProperty(
        "augsys_qn_iter_total",
        ("Total number of augmented system solve iterations used on the "
        "quasi-normal step"))
    augsys_pg_iter_total = Optizelle.createNatProperty(
        "augsys_pg_iter_total",
        ("Total number of augmented system solve iterations used when "
        "projecting the gradient prior to the tangential subproblem"))
    augsys_proj_iter_total = Optizelle.createNatProperty(
        "augsys_proj_iter_total",
        ("Total number of augmented system solve iterations used in the "
        "nullspace projection inside the tangential subproblem"))
    augsys_tang_iter_total = Optizelle.createNatProperty(
        "augsys_tang_iter_total",
        ("Total number of augmented system solve iterations used in the "
        "tangential step"))
    augsys_lmh_iter_total = Optizelle.createNatProperty(
        "augsys_lmh_iter_total",
        ("Total number of augmented system solve iterations used in the "
        "equality multiplier solve"))
    augsys_qn_err = Optizelle.createFloatProperty(
        "augsys_qn_err",
        ("Error in the augmented system solve used on the "
        "quasi-normal step"))
    augsys_pg_err = Optizelle.createFloatProperty(
        "augsys_pg_err",
        ("Error in the augmented system solve used when "
        "projecting the gradient prior to the tangential subproblem"))
    augsys_proj_err = Optizelle.createFloatProperty(
        "augsys_proj_err",
        ("Error in the augmented system solve used in the "
        "nullspace projection inside the tangential subproblem"))
    augsys_tang_err = Optizelle.createFloatProperty(
        "augsys_tang_err",
        ("Error in the augmented system solve used in the "
        "tangential step"))
    augsys_lmh_err = Optizelle.createFloatProperty(
        "augsys_lmh_err",
        ("Error in the augmented system solve used in the "
        "equality multiplier solve"))
    augsys_qn_err_target = Optizelle.createFloatProperty(
        "augsys_qn_err_target",
        ("Target error in the augmented system solve used on the "
        "quasi-normal step"))
    augsys_pg_err_target = Optizelle.createFloatProperty(
        "augsys_pg_err_target",
        ("Target error in the augmented system solve used when "
        "projecting the gradient prior to the tangential subproblem"))
    augsys_proj_err_target = Optizelle.createFloatProperty(
        "augsys_proj_err_target",
        ("Target error in the augmented system solve used in the "
        "nullspace projection inside the tangential subproblem"))
    augsys_tang_err_target = Optizelle.createFloatProperty(
        "augsys_tang_err_target",
        ("Target error in the augmented system solve used in the "
        "tangential step"))
    augsys_lmh_err_target = Optizelle.createFloatProperty(
        "augsys_lmh_err_target",
        ("Target error in the augmented system solve used in the "
        "equality multiplier solve"))
    augsys_iter_total = Optizelle.createNatProperty(
        "augsys_iter_total",
        ("Total number of augmented system solve iterations used in all "
        "solves"))
    augsys_qn_failed = Optizelle.createNatProperty(
        "augsys_qn_failed",
        ("Number of failed quasinormal augmented system solves"))
    augsys_pg_failed = Optizelle.createNatProperty(
        "augsys_pg_failed",
        ("Number of failed projected gradient augmented system solves"))
    augsys_proj_failed = Optizelle.createNatProperty(
        "augsys_proj_failed",
        ("Number of failed nullspace projection augmented system solves"))
    augsys_tang_failed = Optizelle.createNatProperty(
        "augsys_tang_failed",
        ("Number of failed tangential step augmented system solves"))
    augsys_lmh_failed = Optizelle.createNatProperty(
        "augsys_lmh_failed",
        ("Number of failed equality multiplier augmented system solves"))
    augsys_failed_total = Optizelle.createNatProperty(
        "augsys_failed_total",
        ("Total number of failed augmented system solves"))
    g_x = Optizelle.createVectorProperty(
        "g_x",
        ("Equality constraint evaluated at x.  This is used in the quasinormal "
        "step as well as in the computation of the linear Taylor series at x "
        "in the direciton dx_n."))
    norm_gxtyp = Optizelle.createFloatProperty(
        "norm_gxtyp",
        ("A typical norm for norm_gx.  Generally, we just take the value at "
        "the first iteration."))
    gpxdxn_p_gx = Optizelle.createVectorProperty(
        "gpxdxn_p_gx", 
        ("Linear Taylor series at x in the direction dx_n.  This is used both "
        "in the predicted reduction as well as the residual predicted "
        "reduction."))
    gpxdxt = Optizelle.createVectorProperty(
        "gpxdxt",
        ("Derivative of the constraint applied to the tangential step this is "
        "used in the residual predicted reduction."))
    norm_gpxdxnpgx = Optizelle.createVectorProperty(
        "norm_gpxdxnpgx",
        ("Norm of gpxdxn_p_gx.  This is used in the penalty parameter "
        "computation and predicted reduction."))
    dx_n = Optizelle.createVectorProperty(
        "dx_n",
        "Normal step")
    dx_ncp = Optizelle.createVectorProperty(
        "dx_ncp",
        "Cauchy point for normal step")
    dx_t = Optizelle.createVectorProperty(
        "dx_t",
        "(Corrected) tangential step")
    dx_t_uncorrected = Optizelle.createVectorProperty(
        "dx_t_uncorrected",
        "Tangential step prior to correction")
    dx_tcp_uncorrected = Optizelle.createVectorProperty(
        "dx_tcp_uncorrected",
        "Cauchy point for tangential step prior to correction")
    H_dxn = Optizelle.createVectorProperty(
        "H_dxn",
        ("Hessian applied to the normal step.  This is required by W_gradpHdxn "
        "as well as the predicted reduction."))
    W_gradpHdxn = Optizelle.createVectorProperty(
        "W_gradpHdxn",
        ("Quantity grad f(x) + g'(x)*y + H dx_n projected into the null-space ",
        "of the constraints.  This is required in the tangential subproblem "
        "and the predicted reduction."))
    H_dxtuncorrected = Optizelle.createVectorProperty(
        "H_dxtuncorrected",
        ("Hessian applied to the uncorrected tangential step.  This is needed "
        "in the predicted reduction."))
    g_diag = Optizelle.createEnumProperty(
        "g_diag",
        Optizelle.FunctionDiagnostics,
        "Function diagnostics on g")
    y_diag = Optizelle.createEnumProperty(
        "y_diag",
        Optizelle.VectorSpaceDiagnostics,
        "Vector space diagnostics on Y")
    qn_stop = Optizelle.createEnumProperty(
        "qn_stop",
        Optizelle.QuasinormalStop,
        "Reason why the quasinormal problem exited")

def checkT(name,value):
    """Check that we have a state"""
    if not issubclass(type(value),t):
        raise TypeError(
            "The %s argument must have type EqualityConstrained.State.t."
            % (name))
