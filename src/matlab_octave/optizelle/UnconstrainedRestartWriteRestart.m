% Writes a json restart file 
function UnconstrainedRestartWriteRestart(X,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkString('fname',fname);
    checkUnconstrainedStateT('state',state);

    % Write the restart file 
    UnconstrainedRestartWriteRestart_(X,fname,state);
