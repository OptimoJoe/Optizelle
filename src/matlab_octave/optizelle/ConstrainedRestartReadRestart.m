% Reads a json restart file
function state = ConstrainedRestartReadRestart(X,Y,Z,fname,x,y,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkString('fname',fname);

    % Read the restart file
    state = ConstrainedRestartReadRestart_(X,Y,Z,fname,x,y,z);
