% Creates an inequality constrained state 
function self=InequalityConstrainedStateCreate(X,Z,msg,x,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);

    % Create the state
    self=InequalityConstrainedStateCreate_(X,Z,msg,x,z);
