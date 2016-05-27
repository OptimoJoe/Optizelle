% Solves an equality constrained optimization problem
function self=EqualityConstrainedAlgorithmsGetMin(varargin)
    % Check the number of arguments
    if nargin~=5 && nargin~=6
        error(['The getMin function requires either 5 or 6 arguments, ', ...
            sprintf('but %d given.',nargin)]); 
    end

    % Extract the arguments
    X=varargin{1};
    Y=varargin{2};
    msg=varargin{3};
    fns = varargin{4};
    state = varargin{5};
    if nargin==5
        smanip = struct('eval',@(fns,state,loc)state);
    else
        smanip = varargin{6};
    end

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkMessaging('msg',msg);
    checkEqualityConstrainedFunctionsT('fns',fns);
    checkEqualityConstrainedStateT('state',state);
    checkStateManipulator('smanip',smanip);

    % Solve the optimization problem 
    self=EqualityConstrainedAlgorithmsGetMin_(X,Y,msg,fns,state,smanip);
