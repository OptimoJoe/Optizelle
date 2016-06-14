% Release the state in a optimization problem 
function [xs,ys,zs,reals,nats,params]=ConstrainedRestartRelease( ...
    X,Y,Z,state)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkConstrainedStateT('state',state);

    % Do the state release 
    [xs,ys,zs,reals,nats,params]=ConstrainedRestartRelease_( ...
        X,Y,Z,state);
