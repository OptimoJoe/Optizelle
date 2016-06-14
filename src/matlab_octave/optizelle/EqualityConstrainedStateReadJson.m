% Reads equality constrained state parameters from file 
function self=EqualityConstrainedStateReadJson(X,Y,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkString('fname',fname);
    checkEqualityConstrainedStateT('state',state);

    % Read the json file 
    self=EqualityConstrainedStateReadJson_(X,Y,fname,state);
