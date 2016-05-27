% Solves an unconstrained optimization problem
function self=UnconstrainedAlgorithmsGetMin(varargin)
    % Check the number of arguments
    if nargin~=4 && nargin~=5
        error(['The getMin function requires either 4 or 5 arguments, ', ...
            sprintf('but %d given.',nargin)]); 
    end

    % Extract the arguments
    X=varargin{1};
    msg=varargin{2};
    fns = varargin{3};
    state = varargin{4};
    if nargin==4
        smanip = struct('eval',@(fns,state,loc)state);
    else
        smanip = varargin{5};
    end

    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkUnconstrainedFunctionsT('fns',fns);
    checkUnconstrainedStateT('state',state);
    checkStateManipulator('smanip',smanip);

    % Solve the optimization problem 
    self=UnconstrainedAlgorithmsGetMin_(X,msg,fns,state,smanip);
