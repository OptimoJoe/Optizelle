% Optizelle optimization library
%
% Usage:
% global Optizelle;
% setupOptizelle();
%
function setupOptizelle()

global Optizelle

% Intialization is expensive.  Don't do it more than once.
if isempty(Optizelle) 

% Add the optizelle directory to the path, which contains a number of helper
% functions.  In theory, I could just use pwd, but I really want to know where
% *this* file is and not the current directory.
dir = mfilename('fullpath');
dir = dir(1:end-15);
optizelle_dir = sprintf('%s/optizelle',dir);
addpath(optizelle_dir);

% Add jsonlab if we install it in this path
jsonlab_dir = sprintf('%s/jsonlab',dir);
if exist(jsonlab_dir,'dir')
    addpath(jsonlab_dir)
end

% Creates an enumerated type from a cell list of names
createEnum = @(x) cell2struct( ...
    [num2cell(1:length(x)) {@(i)x{i}}], ...
    [x,{'to_string'}],2);

% Reasons we stop truncated CG 
Optizelle.TruncatedStop = createEnum( { ...
    'NotConverged', ...
    'NegativeCurvature', ...
    'RelativeErrorSmall', ...
    'MaxItersExceeded', ...
    'TrustRegionViolated', ...
    'NanOperator', ...
    'NanPreconditioner', ...
    'NonProjector', ...
    'NonSymmetric', ...
    'LossOfOrthogonality', ...
    'InvalidTrustRegionOffset', ...
    'TooManyFailedSafeguard', ...
    'ObjectiveIncrease'} );

% Which algorithm Optizelle.do we use
Optizelle.AlgorithmClass = createEnum( { ...
    'TrustRegion', ...
    'LineSearch', ...
    'UserDefined' } );

% Reasons why we stop the algorithm
Optizelle.OptimizationStop = createEnum( { ...
    'NotConverged', ...
    'GradientSmall', ...
    'StepSmall', ...
    'MaxItersExceeded', ...
    'InteriorPointInstability', ...
    'GlobalizationFailure', ...
    'UserDefined' } );

% Various operators for both Hessian approximations and preconditioners
Optizelle.Operators = createEnum( { ...
    'Identity', ...
    'ScaledIdentity', ...
    'BFGS', ...
    'InvBFGS', ...
    'SR1', ...
    'InvSR1', ...
    'UserDefined' } );
    
% Different kinds of search directions
Optizelle.LineSearchDirection = createEnum( { ...
    'SteepestDescent', ...
    'FletcherReeves', ...
    'PolakRibiere', ...
    'HestenesStiefel', ...
    'BFGS', ...
    'NewtonCG' } );
   
% Different sorts of line searches
Optizelle.LineSearchKind = createEnum( { ...
    'GoldenSection', ...
    'BackTracking', ...
    'TwoPointA', ...
    'TwoPointB' } );
    
% Different points in the optimization algorithm
Optizelle.OptimizationLocation = createEnum( { ...
    'BeginningOfOptimization', ...
    'BeforeInitialFuncAndGrad', ...
    'AfterInitialFuncAndGrad', ...
    'BeforeOptimizationLoop', ...
    'BeginningOfOptimizationLoop', ...
    'BeforeSaveOld', ...
    'BeforeStep', ...
    'BeforeGetStep', ...
    'GetStep', ...
    'AfterStepBeforeGradient', ...
    'AfterGradient', ...
    'BeforeQuasi', ...
    'AfterQuasi', ...
    'AfterCheckStop', ...
    'EndOfOptimizationIteration', ...
    'BeforeLineSearch', ...
    'AfterRejectedTrustRegion', ...
    'AfterRejectedLineSearch', ...
    'BeforeActualVersusPredicted', ...
    'EndOfOptimization' } );

% Different problem classes
Optizelle.ProblemClass = createEnum( { ...
    'Unconstrained', ...
    'EqualityConstrained', ...
    'InequalityConstrained', ...
    'Constrained' } );
    
% Different schemes for adjusting the interior point centrality
Optizelle.CentralityStrategy = createEnum( { ...
    'Constant', ...
    'StairStep', ...
    'PredictorCorrector' } );

% Different function diagnostics on the optimization functions 
Optizelle.FunctionDiagnostics = createEnum( { ...
    'NoDiagnostics', ...
    'FirstOrder', ...
    'SecondOrder' } );

% Different diagnostics on the vector space algebras 
Optizelle.VectorSpaceDiagnostics = createEnum( { ...
    'NoDiagnostics', ...
    'Basic', ...
    'EuclideanJordan' } );

% When and how often we compute our intrusive diagnostics
Optizelle.DiagnosticScheme = createEnum( { ...
    'Never', ...
    'DiagnosticsOnly', ...
    'EveryIteration' } );

% Different kinds of stopping tolerances 
Optizelle.ToleranceKind = createEnum( { ...
    'Absolute', ...
    'Relative'});

% Reasons why the quasinormal problem exited
Optizelle.QuasinormalStop = createEnum( { ...
    'Newton', ...
    'CauchyTrustRegion', ...
    'CauchySafeguard', ...
    'DoglegTrustRegion', ...
    'DoglegSafeguard', ...
    'NewtonTrustRegion', ...
    'NewtonSafeguard', ...
    'Skipped', ...
    'CauchySolved'});

% Different cones used in SQL problems 
Optizelle.Cone = createEnum( { ...
    'Linear', ...
    'Quadratic', ...
    'Semidefinite' } );

%---ScalarValuedFunction0---
% A simple scalar valued function interface, f : X -> R
err_svf=@(x)error(sprintf( ...
    'The %s function is not defined in a ScalarValuedFunction.',x));
Optizelle.ScalarValuedFunction = struct( ...
    'eval',@(x)err_svf('eval'), ...
    'grad',@(x)err_svf('grad'), ...
    'hess_vec',@(x,dx)err_svf('hess_vec'));
%---ScalarValuedFunction1---

%---VectorValuedFunction0---
% A vector valued function interface, f : X -> Y
err_vvf=@(x)error(sprintf( ...
    'The %s function is not defined in a VectorValuedFunction.',x));
Optizelle.VectorValuedFunction = struct( ...
    'eval',@(x)err_vvf('eval'), ...
    'p',@(x,dx)err_vvf('p'), ...
    'ps',@(x,dy)err_vvf('ps'), ...
    'pps',@(x,dx,dy)err_vvf('pps'));
%---VectorValuedFunction1---

%---Operator0---
% A linear operator specification, A : X->Y 
err_op=@(x)error(sprintf( ...
    'The %s function is not defined in an Operator.',x));
Optizelle.Operator = struct( ...
    'eval',@(state,x)err_op('eval'));
%---Operator1---

%---Messaging0---
% Defines how we output messages to the user
Optizelle.Messaging = struct( ...
    'print',@(x)fprintf('%s\n',x), ...
    'error',@(x)fprintf('%s\n',x));
%---Messaging1---

%---StateManipulator0---
% A function that has free reign to manipulate or analyze the state.
Optizelle.StateManipulator = struct('eval',@(fns,state,loc)state);
%---StateManipulator1---

% Vector space for the nonnegative orthant.  For basic vectors in R^m, use this.
Optizelle.Rm = struct( ...
    'init',@(x)x, ...
    'copy',@(x)x, ...
    'scal',@(alpha,x)alpha*x, ...
    'zero',@(x)zeros(size(x)), ...
    'axpy',@(alpha,x,y)alpha*x+y, ...
    'innr',@(x,y)x'*y, ...
    'rand',@(x)randn(size(x)), ...
    'prod',@(x,y)x.*y, ...
    'id',@(x)ones(size(x)), ...
    'linv',@(x,y)y./x, ...
    'barr',@(x)sum(log(x)), ...
    'srch',@(x,y)rm_srch(x,y), ...
    'symm',@(x)x);

% A vector spaces consisting of a finite product of semidefinite,
% quadratic, and linear cones.  This uses the nonsymmetric product
% for the SDP blocks where x o y = xy.  This is not a true Euclidean-Jordan
% algebra, but is sufficient for our purposes.
Optizelle.SQL = struct( ...
    'create',@(types,sizes)sql_create(types,sizes), ...
    'init',@(x)x, ...
    'copy',@(x)x, ...
    'scal',@(alpha,x)sql_scal(alpha,x), ...
    'zero',@(x)sql_zero(x), ...
    'axpy',@(alpha,x,y)sql_axpy(alpha,x,y), ...
    'innr',@(x,y)sql_innr(x,y), ...
    'rand',@(x)sql_rand(x), ...
    'prod',@(x,y)sql_prod(x,y), ...
    'id',@(x)sql_id(x), ...
    'linv',@(x,y)sql_linv(x,y), ...
    'barr',@(x)sql_barr(x), ...
    'srch',@(x,y)sql_srch(x,y), ...
    'symm',@(x)sql_symm(x));

% Converts a vector to a JSON formatted string
Optizelle.json.Serialization.serialize = @serialize;

% Serializes a matlab column vector for the vector space Optizelle.Rm
Optizelle.json.Serialization.serialize( ...
    'register', ...
    @(x,name,iter)strrep(mat2str(x'),' ',', '), ...
    @(x)isvector(x) && isnumeric(x));

% For the moment, do an empty serialization of the SQL vector space 
Optizelle.json.Serialization.serialize( ...
    'register', ...
    @(x,name,iter)savejson('',x), ...
    @(x)isfield(x,'data') && isfield(x,'types') && isfield(x,'sizes'));

% Converts a JSON formatted string to a vector 
Optizelle.json.Serialization.deserialize = @deserialize;

% Deserializes a matlab column vector for the vector space Optizelle.Rm
Optizelle.json.Serialization.deserialize( ...
    'register', ...
    @(x,x_json)str2num(x_json)', ...
    @(x)isvector(x) && isnumeric(x));

% For the moment, do an empty deserialization of the SQL vector space 
Optizelle.json.Serialization.deserialize( ...
    'register', ...
    @(x,x_json)loadjson(x_json), ...
    @(x)isfield(x,'data') && isfield(x,'types') && isfield(x,'sizes'));

%Creates an unconstrained state
Optizelle.Unconstrained.State.t = @UnconstrainedStateCreate;

% All the functions required by an optimization algorithm
Optizelle.Unconstrained.Functions.t= struct( ...
    'f',Optizelle.ScalarValuedFunction, ...
    'PH',Optizelle.Operator);

% Solves an unconstrained optimization problem
Optizelle.Unconstrained.Algorithms.getMin = @UnconstrainedAlgorithmsGetMin;

% Holds restart information
Optizelle.Unconstrained.Restart.X_Vectors = {};
Optizelle.Unconstrained.Restart.Reals = {};
Optizelle.Unconstrained.Restart.Naturals = {};
Optizelle.Unconstrained.Restart.Params = {};

% Release the state in an unconstrained optimization problem 
Optizelle.Unconstrained.Restart.release = @UnconstrainedRestartRelease;

% Capture the state in an unconstrained optimization problem 
Optizelle.Unconstrained.Restart.capture = @UnconstrainedRestartCapture;

% Reads unconstrained state parameters from file 
Optizelle.json.Unconstrained.read = @UnconstrainedStateReadJson;

% Writes a json restart file
Optizelle.json.Unconstrained.write_restart = @UnconstrainedRestartWriteRestart;

% Reads a json restart file
Optizelle.json.Unconstrained.read_restart = @UnconstrainedRestartReadRestart;

%Creates an equality constrained state
Optizelle.EqualityConstrained.State.t = @EqualityConstrainedStateCreate;

% All the functions required by an optimization algorithm
Optizelle.EqualityConstrained.Functions.t= mergeStruct( ...
    Optizelle.Unconstrained.Functions.t, ...
    struct( ...
        'g',Optizelle.VectorValuedFunction, ...
        'PSchur_left',Optizelle.Operator, ...
        'PSchur_right',Optizelle.Operator));

% Solves an equality constrained optimization problem
Optizelle.EqualityConstrained.Algorithms.getMin = ...
    @EqualityConstrainedAlgorithmsGetMin;

% Holds restart information
Optizelle.EqualityConstrained.Restart.X_Vectors = {};
Optizelle.EqualityConstrained.Restart.Y_Vectors = {};
Optizelle.EqualityConstrained.Restart.Reals = {};
Optizelle.EqualityConstrained.Restart.Naturals = {};
Optizelle.EqualityConstrained.Restart.Params = {};

% Release the state in an equality constrained optimization problem 
Optizelle.EqualityConstrained.Restart.release = ...
    @EqualityConstrainedRestartRelease;

% Capture the state in an equality constrained optimization problem 
Optizelle.EqualityConstrained.Restart.capture = ...
    @EqualityConstrainedRestartCapture;

% Reads equality constrained state parameters from file 
Optizelle.json.EqualityConstrained.read = @EqualityConstrainedStateReadJson;

% Writes a json restart file
Optizelle.json.EqualityConstrained.write_restart = ...
    @EqualityConstrainedRestartWriteRestart;

% Reads a json restart file
Optizelle.json.EqualityConstrained.read_restart = ...
    @EqualityConstrainedRestartReadRestart;

%Creates an inequality constrained state
Optizelle.InequalityConstrained.State.t = @InequalityConstrainedStateCreate;

% All the functions required by an optimization algorithm
Optizelle.InequalityConstrained.Functions.t= mergeStruct( ...
    Optizelle.Unconstrained.Functions.t, ...
    struct( ...
        'h',Optizelle.VectorValuedFunction));

% Solves an inequality constrained optimization problem
Optizelle.InequalityConstrained.Algorithms.getMin = ...
    @InequalityConstrainedAlgorithmsGetMin;

% Holds restart information
Optizelle.InequalityConstrained.Restart.X_Vectors = {};
Optizelle.InequalityConstrained.Restart.Z_Vectors = {};
Optizelle.InequalityConstrained.Restart.Reals = {};
Optizelle.InequalityConstrained.Restart.Naturals = {};
Optizelle.InequalityConstrained.Restart.Params = {};

% Release the state in an inequality constrained optimization problem 
Optizelle.InequalityConstrained.Restart.release = ...
    @InequalityConstrainedRestartRelease;

% Capture the state in an inequality constrained optimization problem 
Optizelle.InequalityConstrained.Restart.capture = ...
    @InequalityConstrainedRestartCapture;

% Reads inequality constrained state parameters from file 
Optizelle.json.InequalityConstrained.read = @InequalityConstrainedStateReadJson;

% Writes a json restart file
Optizelle.json.InequalityConstrained.write_restart = ...
    @InequalityConstrainedRestartWriteRestart;

% Reads a json restart file
Optizelle.json.InequalityConstrained.read_restart = ...
    @InequalityConstrainedRestartReadRestart;

%Creates a constrained state
Optizelle.Constrained.State.t = @ConstrainedStateCreate;

% All the functions required by an optimization algorithm
Optizelle.Constrained.Functions.t= ...
    mergeStruct( ...
        mergeStruct( ...
            Optizelle.Unconstrained.Functions.t, ...
            Optizelle.EqualityConstrained.Functions.t), ...
        Optizelle.InequalityConstrained.Functions.t);

% Solves a constrained optimization problem
Optizelle.Constrained.Algorithms.getMin = @ConstrainedAlgorithmsGetMin;

% Holds restart information
Optizelle.Constrained.Restart.X_Vectors = {};
Optizelle.Constrained.Restart.Y_Vectors = {};
Optizelle.Constrained.Restart.Z_Vectors = {};
Optizelle.Constrained.Restart.Reals = {};
Optizelle.Constrained.Restart.Naturals = {};
Optizelle.Constrained.Restart.Params = {};

% Release the state in a constrained optimization problem 
Optizelle.Constrained.Restart.release = @ConstrainedRestartRelease;

% Capture the state in a constrained optimization problem 
Optizelle.Constrained.Restart.capture = @ConstrainedRestartCapture;

% Reads constrained state parameters from file 
Optizelle.json.Constrained.read = @ConstrainedStateReadJson;

% Writes a json restart file
Optizelle.json.Constrained.write_restart = @ConstrainedRestartWriteRestart;

% Reads a json restart file
Optizelle.json.Constrained.read_restart = @ConstrainedRestartReadRestart;
end

end

% Rm line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
% where y > 0
function alpha = rm_srch(x,y)
    % Grab the indices of x less than zero
    i = find(x<0);

    % If we're not moving toward the boundary of the cone in any direction, exit
    if length(i) == 0
        alpha = inf;
        return
    end

    % Find -y / x for the indices where x < 0.  The smallest such number is how
    % far we can travel.
    alpha = min(-y(i)./x(i));
end

% Merges two structures
function ret=mergeStruct(s1,s2)
    % Find the field names and data 
    data = [struct2cell(s1)' struct2cell(s2)'];
    fields = [fieldnames(s1)' fieldnames(s2)'];

    % Figure out the unique elements
    [fields i] = unique(fields);

    % Merge the unique elements
    ret=cell2struct(data(i),fields,2);
end

% Creates an SQL vector
function x = sql_create(types,sizes)
    % Load in Optizelle
    global Optizelle

    % Make sure we have 1-dimensional objects for types and sizes
    if(~isvector(types))
        error('Argument types must be a vector');
    end
    if(~isvector(sizes))
        error('Argument sizes must be a vector');
    end

    % Make sure that types and sizes are the same size
    if length(types) ~= length(sizes)
        error('Arguments types and sizes must be the same length');
    end
                
    % Make sure we have at least one cone.
    if length(types) == 0
        error('A SQL vector requires at least one cone');
    end

    % Make sure that we have something with valid types
    for i=1:length(types)
        if types(i)~=Optizelle.Cone.Linear && ...
           types(i)~=Optizelle.Cone.Quadratic && ...
           types(i)~=Optizelle.Cone.Semidefinite
           error('Parameter types (the first argument) uses invalid cones');
        end
    end

    % Assign the types and sizes
    x.types = types;
    x.sizes = sizes;

    % Create the data 
    for i=1:length(x.sizes)
        if  x.types(i)==Optizelle.Cone.Linear || ...
            x.types(i)==Optizelle.Cone.Quadratic
            x.data{i}=zeros(x.sizes(i),1);
        else
            x.data{i}=zeros(x.sizes(i));
        end
    end
end

% SQL scalar multiply, x <- alpha * x
function x = sql_scal(alpha,x)
    % Loop over all the blocks
    for blk=1:length(x.sizes)
        x.data{blk}=alpha*x.data{blk};
    end
end

% SQL zero, x <- 0
function x = sql_zero(x)
    % Loop over all the blocks
    for blk=1:length(x.sizes)
        x.data{blk}=zeros(size(x.data{blk}));
    end
end

% SQL addition,  y <- alpha * x + y
function y = sql_axpy(alpha,x,y)
    % Loop over all the blocks
    for blk=1:length(x.sizes)
        y.data{blk}=alpha*x.data{blk}+y.data{blk};
    end
end

% SQL inner product, innr <- <x,y>
function z = sql_innr(x,y)
    % Accumulate into z
    z = 0.;

    % Loop over all the blocks
    for blk=1:length(x.sizes)
        z = z + x.data{blk}(:)'*y.data{blk}(:);
    end
end

% SQL random, x <- 0
function x = sql_rand(x)
    % Loop over all the blocks
    for blk=1:length(x.sizes)
        x.data{blk}=randn(size(x.data{blk}));
    end

    % Symmetrize the result
    x = sql_symm(x);
end

% SQL Jordan product, z <- x o y
function z = sql_prod(x,y)
    % Load in Optizelle
    global Optizelle

    % Create a new SQL vector 
    z = sql_create(x.types,x.sizes);

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % z = diag(x) y
        if x.types(blk)==Optizelle.Cone.Linear
            z.data{blk}=x.data{blk}.*y.data{blk};

        % z = [x'y ; x0 ybar + y0 xbar].
        elseif x.types(blk)==Optizelle.Cone.Quadratic
            z.data{blk}=[ ...
                x.data{blk}'*y.data{blk}; ...
                x.data{blk}(1)*y.data{blk}(2:end) + ...
                    y.data{blk}(1)*x.data{blk}(2:end) ];

        % z = xy 
        else 
            z.data{blk}=x.data{blk}*y.data{blk};
        end
    end
end

% SQL identity element, x <- e such that x o e = x
function x = sql_id(x)
    % Load in Optizelle
    global Optizelle

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % z = [1;...;1] 
        if x.types(blk)==Optizelle.Cone.Linear
            x.data{blk}=ones(x.sizes(blk),1);

        % z = [1;0;...;0]
        elseif x.types(blk)==Optizelle.Cone.Quadratic
            x.data{blk}=zeros(x.sizes(blk),1);
            x.data{blk}(1)=1;

        % z = I
        else 
            x.data{blk}=eye(x.sizes(blk));
        end
    end
end

% This applies the inverse of the Schur complement of the Arw operator to a
% vector, inv(Schur(Arw(x)))y.  Note, ybar has size one less than x.
function z = invSchur(x,ybar)
    % z<- 1/x0 ybar + <xbar,ybar> / (x0 ( x0^2 - <xbar,xbar> )) xbar
    z = (1/x(1)) * ybar + ...
        (x(2:end)'*ybar)/(x(1)*(x(1)^2-x(2:end)'*x(2:end))) * x(2:end);
end
        
% SQL Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
function z = sql_linv(x,y)
    % Load in Optizelle
    global Optizelle

    % Create a new SQL vector 
    z = sql_create(x.types,x.sizes);

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % z = inv(Diag(x)) y 
        if x.types(blk)==Optizelle.Cone.Linear
            z.data{blk}=y.data{blk}./x.data{blk};

        % z0 <- y0/(x0-(1/x0) <x_bar,x_bar>) - (1/x0) <x_bar,invSchur(x)(y_bar)>
        % z_bar <- (-y0/x0) invSchur(x)(x_bar) + invSchur(x)(y_bar)
        elseif x.types(blk)==Optizelle.Cone.Quadratic
            invSchur_x_ybar = invSchur(x.data{blk},y.data{blk}(2:end));
            invSchur_x_xbar = invSchur(x.data{blk},x.data{blk}(2:end));
            z.data{blk}(1) = ...
                y.data{blk}(1) / ...
                    (x.data{blk}(1) - ...
                        x.data{blk}(2:end)'*x.data{blk}(2:end) ...
                        / x.data{blk}(1)) - ...
                x.data{blk}(2:end)'*invSchur_x_ybar / x.data{blk}(1);
            z.data{blk}(2:end) = ...
                (-y.data{blk}(1)/x.data{blk}(1)) * invSchur_x_xbar + ...
                invSchur_x_ybar;

        % z = inv(x)y 
        else
            z.data{blk}=x.data{blk}\y.data{blk};
        end
    end
end
        
% SQL barrier function, barr <- barr(x) where x o grad barr(x) = e
function z = sql_barr(x)
    % Load in Optizelle
    global Optizelle

    % This accumulates the barrier's value
    z=0.;

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % z += sum_i log(x_i)
        if x.types(blk)==Optizelle.Cone.Linear
            z = z + sum(log(x.data{blk}));
                
        % z += 0.5 * log(x0^2-<xbar,xbar>)
        elseif x.types(blk)==Optizelle.Cone.Quadratic
            z = z + 0.5 * log( ...
                x.data{blk}(1)^2-x.data{blk}(2:end)'*x.data{blk}(2:end)); 

        % z += log(det(x)).  We compute this by noting that
        % log(det(x)) = log(det(u'u)) = log(det(u')det(u))
        %             = log(det(u)^2) = 2 log(det(u))
        else
            [u p] = chol(x.data{blk});
            
            % Make sure to check if p is not positive definite.  If not, throw
            % and inf.
            if p
                z = inf;
            else
                z = z + 2. * log(prod(diag(u)));
            end
        end
    end
end

% SQL line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
% where y > 0
function alpha = sql_srch(x,y)
    % Load in Optizelle
    global Optizelle

    % Line search parameter
    alpha = inf; 

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % Pointwise, alpha_i = -y_i / x_i.  If this number is positive,
        % then we need to restrict how far we travel.
        if x.types(blk)==Optizelle.Cone.Linear
            alpha0 = rm_srch(x.data{blk},y.data{blk});
                
        % We choose the smallest positive number between:
        % -y0/x0, and the roots of alpha^2 a + alpha b + c
        %
        % where
        %
        % a = x0^2 - ||xbar||^2
        % b = 2x0y0 - 2 <xbar,ybar>
        % c = y0^2 - ||ybar||^2 
        % 
        % Technically, if a is zero, the quadratic formula doesn't
        % apply and we use -c/b instead of the roots.  If b is zero
        % and a is zero, then there's no limit to the line search
        elseif x.types(blk)==Optizelle.Cone.Quadratic

            % Now, first we have to insure that the leading coefficient
            % of the second order cone problem remains nonnegative.
            % This number tells us how far we can step before this
            % is not true.
            if x.data{blk}(1) < 0
                alpha0 = -y.data{blk}(1)/x.data{blk}(1);
            else
                alpha0 = inf;
            end

            % Next, assuming that the leading coefficient is fine,
            % figure out how far we can step before we violate the
            % rest of the SOCP constraint.  This involves solving
            % the quadratic equation from above.
            a = x.data{blk}(1)^2-x.data{blk}(2:end)'*x.data{blk}(2:end);
            b = 2.0*x.data{blk}(1)*y.data{blk}(1) - ...
                2.0*x.data{blk}(2:end)'*y.data{blk}(2:end);
            c = y.data{blk}(1)^2-y.data{blk}(2:end)'*y.data{blk}(2:end);
            alpha0 = [alpha0;roots([a,b,c])];

            % Now, determine the step length.
            i = find(alpha0 > 0);
            alpha0 = min(alpha0(i));

        % We need to find the solution of the generalized eigenvalue
        % problem alpha X v + Y v = 0.  Since Y is positive definite,
        % we want to divide by alpha to get a standard form
        % generalized eigenvalue problem X v = (-1/alpha) Y v.  This
        % means that we solve the problem X v = lambda Y v and then
        % set alpha = -1/lambda as long as lambda is negative.  Note,
        % our Krylov method will converge to lambda from the right,
        % which is going to give an upper bound on alpha.  This is
        % not good for our line search, since we want a lower bound.
        % However, since we get an absolute estimate of the error
        % in lambda, we can just back off of it by a small amount.
        else
            % For whatever reason, eigs is being really flaky and will error
            % that we're not symmetric when we're off of symmetry by 1e-16 or
            % so.  Hence, we force our matrices to be symmetric even though
            % they should be at this point.
            symm=@(x)triu(x)+triu(x,1)';
            lambda = eigs(symm(x.data{blk}),symm(y.data{blk}), ...
                1,'sa',struct('tol',1e-2));
            if lambda < 0
                alpha0 = -1/lambda;
            else
                alpha0 = inf;
            end
        end

        % Update alpha if we need to restrict more
        if alpha0 < alpha
            alpha = alpha0;
        end
    end
end

% Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
% operator
function x = sql_symm(x)
    % Load in Optizelle
    global Optizelle

    % Loop over all the blocks
    for blk=1:length(x.sizes)

        % Find the symmetric part of X, (X+X')/2
        if x.types(blk)==Optizelle.Cone.Semidefinite
            x.data{blk} = (x.data{blk}+x.data{blk}')/2.;
        end
    end
end
