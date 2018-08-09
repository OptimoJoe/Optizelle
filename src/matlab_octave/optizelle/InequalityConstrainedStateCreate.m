% Creates an inequality constrained state
function self=InequalityConstrainedStateCreate(X,Z,x,z)
    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);

    % Create the state
    self=InequalityConstrainedStateCreate_(X,Z,x,z);
