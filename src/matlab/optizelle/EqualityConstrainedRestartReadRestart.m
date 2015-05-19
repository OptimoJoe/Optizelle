% Reads a json restart file 
function state = EqualityConstrainedRestartReadRestart(X,Y,msg,fname,x,y)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkMessaging('msg',msg);
    checkString('fname',fname);

    % Read the restart file 
    state = EqualityConstrainedRestartReadRestart_(X,Y,msg,fname,x,y);
