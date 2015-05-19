% Writes a json restart file 
function ConstrainedRestartWriteRestart(X,Y,Z,msg,fname,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkString('fname',fname);
    checkConstrainedStateT('state',state);

    % Write the restart file 
    ConstrainedRestartWriteRestart_(X,Y,Z,msg,fname,state);
