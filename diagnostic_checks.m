% Run a bunch of diagnostic tests on our functions

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Allocate memory for an initial guess and equality multiplier 
x = [1.2;2.3;3.4];
y = zeros(3,1);
z = zeros(3,1);

% Create an optimization state
state=Optizelle.Constrained.State.t( ...
    Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,x,y,z);

% Modify the state so that we just run our diagnostics and exit
state.dscheme = Optizelle.DiagnosticScheme.DiagnosticsOnly;
state.f_diag = Optizelle.FunctionDiagnostics.SecondOrder;
state.g_diag = Optizelle.FunctionDiagnostics.SecondOrder;
state.h_diag = Optizelle.FunctionDiagnostics.SecondOrder;
state.iter_max = 2;

% Create a bundle of functions
fns=Optizelle.Constrained.Functions.t;
fns.f.eval=@(x)zeros(size(x));
fns.f.grad=@(x)zeros(size(x));
fns.f.hessvec=@(x,dx)zeros(size(x));
fns.g=generate_hyperbolic();
fns.h=generate_logistic(2.0,3.0,4.0);

% Even though this looks like we're solving an optimization problem,
% we're actually just going to run our diagnostics and then exit.
Optizelle.Constrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fns,state);
