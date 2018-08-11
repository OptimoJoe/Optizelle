% Reads unconstrained state parameters from file
function self=UnconstrainedStateReadJson(X,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkString('fname',fname);
    checkUnconstrainedStateT('state',state);

    % Read the json file
    self=UnconstrainedStateReadJson_(X,fname,state);
