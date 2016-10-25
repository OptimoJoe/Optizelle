% Release the state in an unconstrained optimization problem 
function [xs,reals,nats,params]=UnconstrainedRestartRelease(X,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkUnconstrainedStateT('state',state);

    % Do the state release 
    [xs,reals,nats,params]=UnconstrainedRestartRelease_(X,state);
