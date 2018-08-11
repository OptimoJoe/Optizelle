% Reads a json restart file
function state = UnconstrainedRestartReadRestart(X,fname,x)
    % Check our arguments
    checkVectorSpace('X',X);
    checkString('fname',fname);

    % Read the restart file
    state = UnconstrainedRestartReadRestart_(X,fname,x);
