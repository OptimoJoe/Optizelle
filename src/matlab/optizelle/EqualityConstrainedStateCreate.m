% Creates an equality constrained state 
function self=EqualityConstrainedStateCreate(X,Y,msg,x,y)
    % Check our arguments
    checkVectorSpace('X',X)
    checkVectorSpace('Y',Y)
    checkMessaging('msg',msg)

    % Create the state
    self=EqualityConstrainedStateCreate_(X,Y,msg,x,y);
