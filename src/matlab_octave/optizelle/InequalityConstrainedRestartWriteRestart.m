% Writes a json restart file 
function InequalityConstrainedRestartWriteRestart(X,Z,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkString('fname',fname);
    checkInequalityConstrainedStateT('state',state);

    % Write the restart file 
    InequalityConstrainedRestartWriteRestart_(X,Z,fname,state);
