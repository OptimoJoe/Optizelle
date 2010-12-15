% In this test, we solve a simple parameter matching problem using a TR method
% where the approximation Hessian is given by Gauss-Newton

% Fix a random seed for debugging purposes
randn('seed',0);

% Create the true solution
m=25; % Size of the control
u=randn(1,m);

% Create a series of random linear operators
n=20; % Size of the linear operators
A={};
for i=1:m
    A={A{:} randn(n)};
end
B=randn(n);

% Create the right hand sides 
p=100; % The number of right hand sides as well as experimental data points 
b={};
for i=1:p
    b={b{:} randn(n,1)};
end

% Find the true solution based on this control
y=zeros(n,p);
for i=1:p
    y(:,i)=h(A,B,b,i,u);
end

% Create data to test the solution operator.  This is simply the true solution
% with some random noise
d=y;

% Create the initial control
u0=randn(1,m);
u0=u+1e-3*randn(1,m);

% Set parameters for the TR method
eta1=.01;
eta2=.9;
delta=1;
eps_cg=1e-2;
eps_g=1e-6;
eps_d=1e-6;
max_iter=5;
whichH=1;
whichMinv=1;

% Solve the problem
[ustar niters]=pe_test(A,B,b,d,u0,eta1,eta2,delta,eps_cg,eps_g,eps_d, ...
    max_iter,whichH,whichMinv);

% Reconstruct the solution based on this control
ystar=zeros(n,p);
for i=1:p
    ystar(:,i)=h(A,B,b,i,ustar);
end

% Output some diagnostic information
fprintf('The number of iterations required for convergence is: %d\n',niters);
fprintf('The relative error in the control is: %e\n',norm(u-ustar)/norm(u));
fprintf('The relative error in the state is: %e\n',norm(y-ystar)/norm(y));
