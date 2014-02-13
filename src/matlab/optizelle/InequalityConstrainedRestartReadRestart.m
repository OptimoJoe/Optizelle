% Reads a json restart file 
function state = InequalityConstrainedRestartReadRestart(X,Z,msg,fname,x,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);

    % Read the restart file 
    state = InequalityConstrainedRestartReadRestart_(X,Z,msg,fname,x,z);
