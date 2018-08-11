# Enumerated types

class EnumeratedType(object):
    """A generic enumerated type"""
    @classmethod
    def to_string(cls,i):
        """Converts the enumerated type into a string"""
        return filter(lambda (name,value):value==i,
            cls.__dict__.items())[0][0]

class TruncatedStop(EnumeratedType):
    """Reasons we stop truncated CG"""
    NotConverged, \
    NegativeCurvature, \
    RelativeErrorSmall, \
    MaxItersExceeded, \
    TrustRegionViolated, \
    NanOperator, \
    NanPreconditioner, \
    NonProjectorPreconditioner, \
    NonSymmetricPreconditioner, \
    NonSymmetricOperator, \
    LossOfOrthogonality, \
    OffsetViolatesTrustRegion, \
    OffsetViolatesSafeguard, \
    TooManyFailedSafeguard, \
    ObjectiveIncrease \
    = range(15)

class AlgorithmClass(EnumeratedType):
    """Which algorithm class do we use"""
    TrustRegion, \
    LineSearch, \
    UserDefined \
     = range(3)

class OptimizationStop(EnumeratedType):
    """Reasons why we stop the algorithm"""
    NotConverged, \
    GradientSmall, \
    StepSmall, \
    MaxItersExceeded, \
    InteriorPointInstability, \
    GlobalizationFailure, \
    UserDefined \
    = range(7)

class Operators(EnumeratedType):
    """Various operators for both Hessian approximations and preconditioners"""
    Identity, \
    Zero, \
    ScaledIdentity, \
    BFGS, \
    InvBFGS, \
    SR1, \
    InvSR1, \
    UserDefined \
    = range(8)

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
    GoldenSection, \
    BackTracking, \
    TwoPointA, \
    TwoPointB \
    = range(4)

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
    AfterCheckStop, \
    EndOfOptimizationIteration, \
    BeforeLineSearch, \
    AfterRejectedTrustRegion, \
    AfterRejectedLineSearch, \
    BeforeActualVersusPredicted, \
    EndOfOptimization \
    = range(20)

class ProblemClass(EnumeratedType):
    """Different problem classes"""
    Unconstrained, \
    EqualityConstrained, \
    InequalityConstrained, \
    Constrained \
    = range(4)

class FunctionDiagnostics(EnumeratedType):
    """Different function diagnostics on the optimization functions"""
    NoDiagnostics, \
    FirstOrder, \
    SecondOrder \
    = range(3)

class VectorSpaceDiagnostics(EnumeratedType):
    """Different diagnostics on the vector space algebra"""
    NoDiagnostics, \
    Basic, \
    EuclideanJordan\
    = range(3)

class DiagnosticScheme(EnumeratedType):
    """When and how often we compute our intrusive diagnostics"""
    Never, \
    DiagnosticsOnly, \
    EveryIteration \
    = range(3)

class ToleranceKind(EnumeratedType):
    """Different kinds of stopping tolerances"""
    Relative, \
    Absolute \
    = range(2)

class QuasinormalStop(EnumeratedType):
    """Reasons why the quasinormal problem exited"""
    Newton, \
    CauchyTrustRegion, \
    CauchySafeguard, \
    DoglegTrustRegion, \
    DoglegSafeguard, \
    NewtonTrustRegion, \
    NewtonSafeguard, \
    Feasible, \
    CauchySolved, \
    LocalMin, \
    NewtonFailed \
    = range(11)

