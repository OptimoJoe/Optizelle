% Reads constrained state parameters from file 
function self=ConstrainedStateReadJson(X,Y,Z,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkConstrainedStateT('state',state);

    % Read the json file 
    self=ConstrainedStateReadJson_(X,Y,Z,msg,fname,state);
