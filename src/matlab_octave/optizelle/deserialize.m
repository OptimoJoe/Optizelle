% Converts a JSON formatted string to a vector
function x_out = deserialize(varargin)
    % Keep track of the serialization functions and checks
    persistent sfns;
    persistent schecks;

    % Do some initialization to make sure we have cell arrays
    if isempty(sfns)
        sfns = {};
        schecks = {};
    end

    % Determine if we're registering a new serialization function
    if nargin==3
        % Grab the arguments
        mode = varargin{1};
        sfn = varargin{2};
        scheck = varargin{3};

        % Check the arguments
        checkString('mode',mode);
        checkFunction('deserialize',sfn);
        checkFunction('check',scheck);

        % Check that we're actually registering
        if ~strcmp(mode,'register')
            error('The only valid mode is register');
        end

        % Register the new function and check
        sfns = {sfns{:},sfn};
        schecks = {schecks{:},scheck};

    % Call the serialization on the vector
    elseif nargin==2
        % Grab the arguments
        x = varargin{1};
        x_json = varargin{2};

        % Try to deserialize the vector
        for i=1:length(sfns)
            if schecks{i}(x)
                x_out = sfns{i}(x,x_json);
                break;
            end
        end

        % Check if we couldn't serialize
        if length(schecks) == 0  ||  ~schecks{i}(x)
            error('Unable to find a suitable deserialization function');
        end

    % Throw an error if we have the wrong number of argumnets
    else
        error(['The serialize function must be called either in ' ...
            'registration mode (three arguments) or deserialize mode ' ...
            '(two arguments)']);
    end
end
