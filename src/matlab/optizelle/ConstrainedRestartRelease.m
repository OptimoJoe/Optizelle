% Release the state in a optimization problem 
function [xs,ys,zs,reals,nats,params]=ConstrainedRestartRelease( ...
    X,Y,Z,msg,state)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkConstrainedStateT('state',state);

    % Do the state release 
    [xs,ys,zs,reals,nats,params]=ConstrainedRestartRelease_( ...
        X,Y,Z,msg,state);
