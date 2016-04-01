% Demonstrates how to cache various calculations when solving an unconstrained
% optimization problem and then verifies that this caching works correctly

% Grab Optizelle 
global Optizelle;
setupOptizelle();

% Set the name of the optimization parameters 
pname = 'reduced_space.json';

% Setup some diagnostics that verify that we're caching properly
global diagnostics
diagnostics.used_cached_objective = 0;
diagnostics.state_factorization_cached = 0;
diagnostics.hessian_cached = 0;
diagnostics.hessian_factorization_cached = 0;

% Grab the parameters
params = generate_params();

% Create a bundle of functions
fns=Optizelle.Unconstrained.Functions.t;
[fns.f fns.PH phi] = generate_reduced_space(params);

% Create an optimization state
x = [ 2.; 2.];
state = Optizelle.Unconstrained.State.t(Optizelle.Rm,Optizelle.Messaging,x);

% Read the parameters from file
state = Optizelle.json.Unconstrained.read( ...
    Optizelle.Rm,Optizelle.Messaging,pname,state);

% Solve the optmization problem
state = Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging,fns,state);

% Find the state solution based on the material we just solved for
u = phi(state.x);

% Plot the final result
if 0
    plot( ...
        params.omega,u, ...
        params.omega,params.d, ...
        params.omega,params.u(params.omega));
    legend('Computed State','Data','Exact')
end

% Validate that we cached correctly
if diagnostics.state_factorization_cached ~= state.iter
    error('Missed a cached factorization of the state equations');
end
if diagnostics.used_cached_objective ~= 1
    error('Missed a cached objective evaluation');
end
if diagnostics.hessian_cached ~= state.iter-1
    error('Missed a cached Hessian');
end
if diagnostics.hessian_factorization_cached ~= state.iter-1
    error('Missed a cached factorization of the Hessian');
end
