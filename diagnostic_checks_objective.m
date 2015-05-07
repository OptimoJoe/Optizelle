% Run a bunch of diagnostic tests on our functions

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Generate a random MLP 
ninput = 3;
nhidden = 20;
nsamples = 5;
alpha = randn(nhidden,1);
A = randn(nhidden,ninput);
b = randn(nhidden,1);

% Generate some random data
x = randn(ninput,nsamples);
y = randn(nsamples,1);

% Allocate memory for an initial guess
xx = randn(nhidden+nhidden*ninput+nhidden,1);

% Create an optimization state
state=Optizelle.Unconstrained.State.t( ...
    Optizelle.Rm,Optizelle.Messaging,xx);

% Modify the state so that we just run our diagnostics and exit
state.dscheme = Optizelle.DiagnosticScheme.DiagnosticsOnly;
state.f_diag = Optizelle.FunctionDiagnostics.SecondOrder;
state.iter_max = 2;

% Grab some lenses
lens=generate_lenses(ninput,nhidden);

% Generate the objective function
fns.f = generate_objective(generate_hyperbolic(),lens,x,y);

% Even though this looks like we're solving an optimization problem,
% we're actually just going to run our diagnostics and then exit.
Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging,fns,state);

% Now, test the parametrization function
fns.f = generate_parametrization(generate_hyperbolic(),lens,x);
Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging,fns,state);
