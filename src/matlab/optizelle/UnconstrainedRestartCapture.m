% Capture the state in an unconstrained optimization problem 
function state=UnconstrainedRestartCapture(X,msg,state,xs,reals,nats,params)
    % Check our arguments
    checkVectorSpace('X',X);
    checkMessaging('msg',msg);
    checkUnconstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture 
    state=UnconstrainedRestartCapture_(X,msg,state,xs,reals,nats,params);
