% This finds the minimum of the Rosenbrock function.

% Make sure that Optizelle is in the path
addpath('../../');

% Create the vector spaces
clear X;
X.copy=@(x)x;
X.scal=@(alpha,x)alpha*x;
X.zero=@(x)zeros(size(x));
X.axpy=@(alpha,x,y)alpha*x+y;
X.innr=@(x,y)x'*y;

% Create a bundle of vector spaces
clear VS;
VS.X = X;

% Create the functions

% Define the Rosenbrock function where
% 
% f(x,y)=(1-x)^2+100(y-x^2)^2
%
clear f;
f.eval=@(x)(1-x(1))^2+100*(x(2)-x(1)^2)^2;
f.grad=@(x) ...
    [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); ...
    200*(x(2)-x(1)^2)];
f.hessvec=@(x,dx) ...
    [(1200*(x(1)^2)-400*x(2)+2)*dx(1)-400*x(1)*dx(2); ...
    -400*x(1)*dx(1) + 200*dx(2)];

% Create the bundle of functions
clear fns;
fns.f=f;

% Set a starting guess
pts.x=[-1.2;0.5];

% Optimize the problem 
sol=Optizelle(VS,fns,pts,'rosenbrock.param');
sol.x
