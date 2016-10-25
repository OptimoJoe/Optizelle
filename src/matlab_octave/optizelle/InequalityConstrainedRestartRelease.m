% Release the state in an inequality optimization problem 
function [xs,zs,reals,nats,params]=InequalityConstrainedRestartRelease( ...
    X,Z,state)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkInequalityConstrainedStateT('state',state);

    % Do the state release 
    [xs,zs,reals,nats,params]=InequalityConstrainedRestartRelease_( ...
        X,Z,state);
