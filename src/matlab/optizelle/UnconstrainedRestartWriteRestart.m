% Writes a json restart file 
function UnconstrainedRestartWriteRestart(X,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkUnconstrainedStateT('state',state);

    % Write the restart file 
    UnconstrainedRestartWriteRestart_(X,msg,fname,state);
