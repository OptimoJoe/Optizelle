% Solves a constrained optimization problem
function self=ConstrainedAlgorithmsGetMin(varargin)
    % Check the number of arguments
    if nargin~=6 && nargin~=7
        error(['The getMin function requires either 5 or 6 arguments, ', ...
            sprintf('but %d given',nargin)]);
    end

    % Extract the arguments
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    msg=varargin{4};
    fns = varargin{5};
    state = varargin{6};
    if nargin==6
        smanip = struct('eval',@(fns,state,loc)state);
    else
        smanip = varargin{7};
    end

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkConstrainedFunctionsT('fns',fns);
    checkConstrainedStateT('state',state);
    checkStateManipulator('smanip',smanip);

    % Solve the optimization problem
    self=ConstrainedAlgorithmsGetMin_(X,Y,Z,msg,fns,state,smanip);
