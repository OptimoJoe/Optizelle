% Creates an unconstrained state
function self=UnconstrainedStateCreate(X,x)
    % Check our arguments
    checkVectorSpace('X',X);

    % Create the state
    self=UnconstrainedStateCreate_(X,x);
