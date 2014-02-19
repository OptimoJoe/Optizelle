% Reads a json restart file 
function state = UnconstrainedRestartReadRestart(X,msg,fname,x)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkString('fname',fname);

    % Read the restart file 
    state = UnconstrainedRestartReadRestart_(X,msg,fname,x);
