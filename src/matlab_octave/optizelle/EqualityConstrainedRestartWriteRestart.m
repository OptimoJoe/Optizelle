% Writes a json restart file
function EqualityConstrainedRestartWriteRestart(X,Y,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkString('fname',fname);
    checkEqualityConstrainedStateT('state',state);

    % Write the restart file
    EqualityConstrainedRestartWriteRestart_(X,Y,fname,state);
