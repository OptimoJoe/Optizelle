% Capture the state in an inequality constrained optimization problem 
function state=InequalityConstrainedRestartCapture( ...
    X,Z,msg,state,xs,zs,reals,nats,params)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkInequalityConstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkVectors('zs',zs);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture 
    state=InequalityConstrainedRestartCapture_( ...
        X,Z,msg,state,xs,zs,reals,nats,params);
