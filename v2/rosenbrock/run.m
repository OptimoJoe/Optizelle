% Add pelab to the path
addpath('..');

% Setup the parameters for the optimization
params=pelab();
params.F='myfunc';
params.G='mygrad';
params.H='myhess';
params.StateManipulator='smanip';
params.iter_max=250;
params.H_type='External';
params.algorithm_class='TrustRegion';
%params.dir='LimitedMemoryBFGS';
%params.stored_history=2;
%params.kind='BackTracking';
params.eps_s=1e-10;
params.eps_g=1e-12;
params.eps_krylov=1e-8;
%params.linesearch_iter_max=20;
%params.alpha=.01;
params.verbose=2;

% Choose an initial guess
x=[-1.2;1];

% Run the optimization
[y why_stop]=pelab(params,x)

