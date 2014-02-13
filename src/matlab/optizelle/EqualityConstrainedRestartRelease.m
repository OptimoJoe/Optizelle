% Release the state in an unconstrained optimization problem 
function [xs,ys,reals,nats,params]=EqualityConstrainedRestartRelease( ...
    X,Y,msg,state)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkMessaging('msg',msg);
    checkEqualityConstrainedStateT('state',state);

    % Do the state release 
    [xs,ys,reals,nats,params]=EqualityConstrainedRestartRelease_(X,Y,msg,state);
