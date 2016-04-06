% Demonstrates how to cache various calculations when solving an equality
% constrained optimization problem and then verifies that this caching works
% correctly

% Grab Optizelle 
global Optizelle;
setupOptizelle();

% Set the name of the optimization parameters 
pname = 'full_space.json';

% Setup some diagnostics that verify that we're caching properly
global diagnostics
diagnostics.first_derivative_cached = 0;
diagnostics.second_derivative_cached = 0;
diagnostics.factorization_cached = 0;

% Grab the parameters
params = generate_params();

% Create a bundle of functions
fns=Optizelle.EqualityConstrained.Functions.t;
[fns.f fns.g fns.PSchur_left phi] = generate_full_space(params);

% Create an optimization state
x = zeros(params.nx+2,1);
    x(params.idx.k) = [ 1.; 1.];
    x(params.idx.u) = phi(x(params.idx.k)); 
y = zeros(params.nx,1);
state = Optizelle.EqualityConstrained.State.t( ...
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,x,y);

% Read the parameters from file
state = Optizelle.json.EqualityConstrained.read( ...
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,pname,state);

% Solve the optmization problem
state = Optizelle.EqualityConstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fns,state);

% Extract the state solution 
u = state.x(params.idx.u); 

% Plot the final result
if 0
    plot( ...
        params.omega,u, ...
        params.omega,params.d, ...
        params.omega,params.u(params.omega));
    legend('Computed State','Data','Exact')
end

% Validate that we cached correctly
if diagnostics.factorization_cached ~= state.glob_iter_total+1
    error('Missed a cached factorization');
end
if diagnostics.first_derivative_cached ~= state.glob_iter_total+1 
    error('Missed a cached first derivative');
end
if diagnostics.second_derivative_cached ~= state.iter-1
    error('Missed a cached second derivative');
end
