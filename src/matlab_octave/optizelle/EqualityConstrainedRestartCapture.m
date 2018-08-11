% Capture the state in an equality constrained optimization problem
function state=EqualityConstrainedRestartCapture( ...
    X,Y,state,xs,ys,reals,nats,params)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkEqualityConstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkVectors('ys',ys);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture
    state=EqualityConstrainedRestartCapture_( ...
        X,Y,state,xs,ys,reals,nats,params);
