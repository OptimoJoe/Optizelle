% Reads a json restart file 
function state = ConstrainedRestartReadRestart(X,Z,msg,fname,x,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);

    % Read the restart file 
    state = ConstrainedRestartReadRestart_(X,Y,Z,msg,fname,x,y,z);
