% Release the state in an unconstrained optimization problem 
function [xs,reals,nats,params]=UnconstrainedRestartRelease(X,msg,state)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkUnconstrainedStateT('state',state);

    % Do the state release 
    [xs,reals,nats,params]=UnconstrainedRestartRelease_(X,msg,state);
