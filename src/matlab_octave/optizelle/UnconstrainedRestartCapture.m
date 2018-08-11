% Capture the state in an unconstrained optimization problem
function state=UnconstrainedRestartCapture(X,state,xs,reals,nats,params)
    % Check our arguments
    checkVectorSpace('X',X);
    checkUnconstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture
    state=UnconstrainedRestartCapture_(X,state,xs,reals,nats,params);
