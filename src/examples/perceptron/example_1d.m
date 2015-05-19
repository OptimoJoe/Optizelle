% Interpolate on a simple 1-D problem 
function example_1d()

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Set the name of the parameter file
fname = 'example_1d.json';

% Set whether we're plotting or not
do_plot = 0;

% Set the size of the problem 
ninput = 1;
nhidden = 10;
nsamples = 20;

% Generate some random data
x = rand(ninput,nsamples);
true_fn = @(x)cos(5*x);
y = true_fn(x); 

% Allocate memory for an initial guess
xx = randn(nhidden+nhidden*ninput+nhidden,1);

% Create an optimization state
state=Optizelle.Unconstrained.State.t( ...
    Optizelle.Rm,Optizelle.Messaging,xx);

% Read the parameters from file
state=Optizelle.json.Unconstrained.read(Optizelle.Rm,Optizelle.Messaging,...
    fname,state);

% Grab some lenses
lens=generate_lenses(ninput,nhidden);

% Generate the objective function
fns=Optizelle.Unconstrained.Functions.t;
fns.f = generate_objective(generate_hyperbolic(),lens,x,y);

% Interpolate our data 
state=Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging,fns,state);

% Generate an interpolatory function based on this
ff = generate_interpolant(generate_hyperbolic(),lens,state.x);

% Return if we're not plotting
if ~do_plot
    return;
end

% Plot the result
x_uniform=0:.01:1;
plot(x_uniform,true_fn(x_uniform),'*',x_uniform,ff.eval(x_uniform),'x');
