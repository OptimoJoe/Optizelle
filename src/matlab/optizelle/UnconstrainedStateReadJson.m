% Reads unconstrained state parameters from file 
function self=UnconstrainedStateReadJson(X,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkUnconstrainedStateT('state',state);

    % Read the json file 
    self=UnconstrainedStateReadJson_(X,msg,fname,state);
