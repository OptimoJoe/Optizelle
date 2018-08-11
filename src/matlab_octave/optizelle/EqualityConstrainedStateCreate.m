% Creates an equality constrained state
function self=EqualityConstrainedStateCreate(X,Y,x,y)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);

    % Create the state
    self=EqualityConstrainedStateCreate_(X,Y,x,y);
