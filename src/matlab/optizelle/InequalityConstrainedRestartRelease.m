% Release the state in an inequality optimization problem 
function [xs,zs,reals,nats,params]=InequalityConstrainedRestartRelease( ...
    X,Z,msg,state)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkInequalityConstrainedStateT('state',state);

    % Do the state release 
    [xs,zs,reals,nats,params]=InequalityConstrainedRestartRelease_( ...
        X,Z,msg,state);
