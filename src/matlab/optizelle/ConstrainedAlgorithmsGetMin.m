% Solves a constrained optimization problem
function self=ConstrainedAlgorithmsGetMin(varargin)
    % Check the number of arguments
    if nargin~=6 && nargin~=7
        error(['The getMin function requires either 5 or 6 arguments, ', ...
            sprintf('but %d given.',nargin)]); 
    end

    % Extract the arguments
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    msg=varargin{4};
    if nargin==6
        fns = varargin{5};
        state = varargin{6};
        smanip = getStateManipulator(); 
    else
        fns = varargin{6};
        state = varargin{7};
        smanip = varargin{5};
    end

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkStateManipulator('smanip',smanip);
    checkConstrainedFunctionsT('fns',fns);
    checkConstrainedStateT('state',state);

    % Solve the optimization problem 
    self=ConstrainedAlgorithmsGetMin_(X,Y,Z,msg,smanip,fns,state);
