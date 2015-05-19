% Converts a vector to a JSON formatted string
function x_json = serialize(varargin)
    % Keep track of the serialization functions and checks
    persistent sfns;
    persistent schecks;

    % Do some initialization to make sure we have cell arrays
    if isempty(sfns)
        sfns = {};
        schecks = {};
    end

    % Determine if we're registering a new serialization function
    if nargin==3 && ischar(varargin{1}) && strcmp(varargin{1},'register')
        % Grab the arguments
        mode = varargin{1};
        sfn = varargin{2};
        scheck = varargin{3};

        % Check the arguments
        checkString('mode',mode);
        checkFunction('serialize',sfn);
        checkFunction('check',scheck);

        % Register the new function and check
        sfns = {sfns{:},sfn};
        schecks = {schecks{:},scheck};

    % Call the serialization on the vector
    elseif nargin==3
        % Grab the arguments
        x = varargin{1};
        name = varargin{2};
        iter = varargin{3};
        
        % Check the arguments
        checkString('name',name);
        checkNatural('iter',iter);

        % Try to serialize the vector
        for i=1:length(sfns)
            if schecks{i}(x)
                x_json = sfns{i}(x,name,iter);
                break;
            end
        end
        
        % Check if we couldn't serialize
        if length(schecks) == 0  ||  ~schecks{i}(x)
            error('Unable to find a suitable serialization function');
        end

    % Throw an error if we have the wrong number of argumnets
    else
        error(['The serialize function must be called either in ' ...
            'registration mode (mode,serialize,check) or serialize mode ' ...
            '(x,name,iter)']);
    end
end
