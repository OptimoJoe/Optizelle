% Tests our ability to throw a native exception and catch it
function objective()
    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Allocate memory for an initial guess
    x = [1.2,2.3];

    % Create an optimization state
    state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x);

    % Create a bundle of functions
    fns=Optizelle.Unconstrained.Functions.t;
    fns.f=Objective();

    % Try to catch the rror
    msg = '';
    try
        state = Optizelle.Unconstrained.Algorithms.getMin( ...
            Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);
    catch e
        msg = e.message;
    end
    if length(msg)==0
        error('Failed to catch an error');
    end
end

% Define a function that just throws an exception
function self = Objective()
    % Evaluation
    self.eval = @myeval;

    % Gradient
    self.grad = @mygrad;

    % Hessian-vector product
    self.hessvec = @myhessvec;
end
function z = myeval(x)
    error('Evaluation');
end
function grad = mygrad(x)
    error('Gradient');
end
function hv = myhessvec(x,dx)
    error('Hessian-vector product')
end
