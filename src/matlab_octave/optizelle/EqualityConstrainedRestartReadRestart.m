% Reads a json restart file 
function state = EqualityConstrainedRestartReadRestart(X,Y,fname,x,y)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkString('fname',fname);

    % Read the restart file 
    state = EqualityConstrainedRestartReadRestart_(X,Y,fname,x,y);
