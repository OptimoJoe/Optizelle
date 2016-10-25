% Writes a json restart file 
function ConstrainedRestartWriteRestart(X,Y,Z,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkString('fname',fname);
    checkConstrainedStateT('state',state);

    % Write the restart file 
    ConstrainedRestartWriteRestart_(X,Y,Z,fname,state);
