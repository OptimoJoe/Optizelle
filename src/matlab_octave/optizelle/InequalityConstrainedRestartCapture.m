% Capture the state in an inequality constrained optimization problem 
function state=InequalityConstrainedRestartCapture( ...
    X,Z,state,xs,zs,reals,nats,params)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkInequalityConstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkVectors('zs',zs);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture 
    state=InequalityConstrainedRestartCapture_( ...
        X,Z,state,xs,zs,reals,nats,params);
