% Reads inequality constrained state parameters from file
function self=InequalityConstrainedStateReadJson(X,Z,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkString('fname',fname);
    checkInequalityConstrainedStateT('state',state);

    % Read the json file
    self=InequalityConstrainedStateReadJson_(X,Z,fname,state);
