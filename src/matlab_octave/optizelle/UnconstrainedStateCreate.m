% Creates an unconstrained state 
function self=UnconstrainedStateCreate(X,msg,x)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);

    % Create the state
    self=UnconstrainedStateCreate_(X,msg,x);
