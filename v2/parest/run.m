% Add in the appropriate optimization routines
addpath('..');

% Setup the problem specifications
randn('seed',0) % Random seed
n=30;		% Size of the state (y)
m=40;		% Size of the control (u)
p=20;		% Number of experiments
beta=0e-1;	% Set the regularization parameter
data_err=0e-5;	% Amount of error in the data
start_err=1e-0; % Amount of error in the starting solution
do_check=1;	% Do a finite difference check?

% Setup the optimization parameters
params=pelab();	
%params.StateManipulator='smanip';
params.iter_max=200;
params.H_type='External';
params.algorithm_class='TrustRegion';
%params.algorithm_class='LineSearch';
%params.dir='LimitedMemoryBFGS';
%params.stored_history=3;
%params.kind='BackTracking';
params.eps_s=1e-10;
params.eps_g=1e-10;
params.eps_krylov=1e-8;
params.krylov_iter_max=200;
%params.linesearch_iter_max=20;
%params.alpha=.01;
params.verbose=1;

% Setup the problem

% Create the norm
mynorm=@(x)sqrt(innr(x,x));

% Generate m random matrices. 
clear A; 
for i=1:m
    A{i}=randn(n);
end

% Generate one extra random matrix for the linear system
B=randn(n);

% Generate the rhs
clear b;
for i=1:p
    b{i}=randn(n,1);
end

% Get the solution operators
h=gen_solution(A,B,b);

% Get a true solution 
u_true=randn(m,1);

% Generate the data
clear d;
for i=1:p
    d{i}=h.f(i,u_true);
    noise=randn(n,1);
    noise=noise/mynorm(d{i})*data_err;
    d{i}=d{i}+noise;
end

% Get the solution operators
r=gen_residual(d);

% Create a list of sample vectors
vecs.u=randn(m,1);
vecs.y=randn(n,1);
vecs.z=randn(n,1);

% Get a starting solution
%u_init=randn(m,1);
u_init=u_true+randn(m,1)/mynorm(u_true)*start_err;

% Setup Tikhonov Regularization
reg_F=@(u).5*beta*mynorm(u-u_init)^2;
reg_G=@(u)beta*(u-u_init);
reg_H=@(u,du)beta*du;

% Build functions for the objective, gradient, and Hessian
F=@(u)getObjective(r,h,vecs,u)+reg_F(u);
G=@(u)getGradient(r,h,vecs,u)+reg_G(u);
%H=@(u,du)getGaussNewton(r,h,vecs,u,du)+reg_H(u,du);
H=@(u,du)getNewton(r,h,vecs,u,du)+reg_H(u,du);
params.F=F;
params.G=G;
params.H=H;

% Do a finite difference test
if do_check 
    eta=randn(m,1);
    gradientCheck(F,G,u_init,eta);
    fpCheck(h,u_init,eta,vecs.y)
    hessianCheck(G,H,u_init,eta)
end

% Solve the optimization problem
[u_sol why_stop]=pelab(params,u_init);

fprintf('The relative error in the parameters is: %e\n', ...
    mynorm(u_true-u_sol)/(1e-16+mynorm(u_true)));

return
