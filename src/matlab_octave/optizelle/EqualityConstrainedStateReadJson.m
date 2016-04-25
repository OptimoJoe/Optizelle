% Reads equality constrained state parameters from file 
function self=EqualityConstrainedStateReadJson(X,Y,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkEqualityConstrainedStateT('state',state);

    % Read the json file 
    self=EqualityConstrainedStateReadJson_(X,Y,msg,fname,state);
