% Creates a constrained state 
function self=ConstrainedStateCreate(X,Y,Z,msg,x,y,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);

    % Create the state
    self=ConstrainedStateCreate_(X,Y,Z,msg,x,y,z);
