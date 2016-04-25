% Writes a json restart file 
function InequalityConstrainedRestartWriteRestart(X,Z,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkInequalityConstrainedStateT('state',state);

    % Write the restart file 
    InequalityConstrainedRestartWriteRestart_(X,Z,msg,fname,state);
