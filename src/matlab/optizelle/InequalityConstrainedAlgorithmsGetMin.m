% Solves an inequality constrained optimization problem
function self=InequalityConstrainedAlgorithmsGetMin(varargin)
    % Check the number of arguments
    if nargin~=5 && nargin~=6
        error(['The getMin function requires either 5 or 6 arguments, ', ...
            sprintf('but %d given.',nargin)]); 
    end

    % Extract the arguments
    X=varargin{1};
    Z=varargin{2};
    msg=varargin{3};
    if nargin==5
        fns = varargin{4};
        state = varargin{5};
        smanip = getStateManipulator(); 
    else
        fns = varargin{5};
        state = varargin{6};
        smanip = varargin{4};
    end

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkStateManipulator('smanip',smanip);
    checkInequalityConstrainedFunctionsT('fns',fns);
    checkInequalityConstrainedStateT('state',state);

    % Solve the optimization problem 
    self=InequalityConstrainedAlgorithmsGetMin_(X,Z,msg,smanip,fns,state);
