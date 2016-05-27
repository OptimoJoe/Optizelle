% Reads inequality constrained state parameters from file 
function self=InequalityConstrainedStateReadJson(X,Z,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkInequalityConstrainedStateT('state',state);

    % Read the json file 
    self=InequalityConstrainedStateReadJson_(X,Z,msg,fname,state);
