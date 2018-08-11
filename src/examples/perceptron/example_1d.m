% Interpolate on a simple 1-D problem

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Set the name of the parameter file
fname = 'example_1d.json';

% Set whether we're plotting or not
do_plot = 0;

% Set the size of the problem
ninput = 1;
nhidden = 3;
nsamples = 20;

% Grab some lenses
[lens idx]=generate_lenses(ninput,nhidden);

% Generate some random data
x = rand(ninput,nsamples);
true_fn = @(x)cos(2*x);
y = true_fn(x);

% Generate scalings based on this data
scaling = generate_scaling(x,y);

% Allocate memory for an initial guess
xx = randn(idx.size,1);

% Create an optimization state
state=Optizelle.Unconstrained.State.t(Optizelle.Rm,xx);

% Read the parameters from file
state=Optizelle.json.Unconstrained.read(Optizelle.Rm,fname,state);

% Generate the objective function
fns=Optizelle.Unconstrained.Functions.t;
fns.f = generate_objective(generate_hyperbolic(),lens,x,y,scaling);

% Interpolate our data
state=Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

% Generate an interpolatory function based on this
ff = generate_interpolant(generate_hyperbolic(),lens,state.x,scaling);

% Return if we're not plotting
if ~do_plot
    return;
end

% Plot the result
x_uniform=0:.01:1;
plot(x_uniform,true_fn(x_uniform),'*',x_uniform,ff.eval(x_uniform),'x');
