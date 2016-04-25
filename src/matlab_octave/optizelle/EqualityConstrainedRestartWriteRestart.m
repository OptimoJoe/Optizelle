% Writes a json restart file 
function EqualityConstrainedRestartWriteRestart(X,Y,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkEqualityConstrainedStateT('state',state);

    % Write the restart file 
    EqualityConstrainedRestartWriteRestart_(X,Y,msg,fname,state);
