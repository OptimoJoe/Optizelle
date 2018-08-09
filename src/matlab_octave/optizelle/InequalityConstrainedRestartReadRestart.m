% Reads a json restart file
function state = InequalityConstrainedRestartReadRestart(X,Z,fname,x,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkString('fname',fname);

    % Read the restart file
    state = InequalityConstrainedRestartReadRestart_(X,Z,fname,x,z);
