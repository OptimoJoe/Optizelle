% Set the random seed
randn('seed',0);

% Set some parameters
m=50; % Size of the state
n=50; % Size of the control
p=100; % Number of experiments
eps_g=1e-6;
eps_d=1e-6;
eps_f=1e-6;
eps_cg=1e-2;
max_cg_iter=5;
eta1=.25;
eta2=.75;
delta=10;
stored_history=10;
history_reset=5;
max_iter=100;

% Create the linear operators
A={};
for i=1:n
    A={A{:} randn(m)};
end
B=randn(m);

% Create the sources 
b={};
for i=1:p
    b={b{:} randn(m,1)};
end

% Create the true solution 
u=randn(n,1);

% Create an initial solution
u0=u+randn(n,1)*1e-1/norm(u,2);

% Create the complete operator based on this control
C=zeros(m);
for i=1:n
    C=C+A{i}*u(i);
end
C=C+B;

% Create the data
d={};
for i=1:p
    d={d{:} C\b{i}+randn(m,1)*1e-4/norm(b{i},2)};
end

% Optimize
u_sol=parest(A,B,b,d,u0,'GaussNewton','Identity',eps_g,eps_d,eps_f,eps_cg, ...
    max_cg_iter,eta1,eta2,delta,stored_history,history_reset,max_iter);

% Check the relative error in the parameters
fprintf('The relative error in the parameters is: %e\n', ...
    norm(u-u_sol)/(1+norm(u)));

% Create an operator based on this solved set of parameters 
C_sol=zeros(m);
for i=1:n
    C_sol=C+A{i}*u_sol(i);
end
C_sol=C_sol+B;
