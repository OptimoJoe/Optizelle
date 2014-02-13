% Capture the state in a constrained optimization problem 
function state=ConstrainedRestartCapture( ...
    X,Y,Z,msg,state,xs,ys,zs,reals,nats,params)

    % Check our arguments
    checkVectorSpace('X',X);
    checkVectorSpace('Y',Y);
    checkVectorSpace('Z',Z);
    checkMessaging('msg',msg);
    checkConstrainedStateT('state',state);
    checkVectors('xs',xs);
    checkVectors('ys',ys);
    checkVectors('zs',zs);
    checkReals('reals',reals);
    checkNaturals('nats',nats);
    checkParams('params',params);

    % Do the state capture 
    state=ConstrainedRestartCapture_( ...
        X,Y,Z,msg,state,xs,ys,zs,reals,nats,params);
