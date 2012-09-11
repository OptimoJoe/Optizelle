% This finds the minimum of the Rosenbrock function with simple inequalities.

% Make sure peopt is in the path
addpath('../../');

% Create the vector spaces
clear X;
X.copy=@(x)x;
X.scal=@(alpha,x)alpha*x;
X.zero=@(x)zeros(size(x));
X.axpy=@(alpha,x,y)alpha*x+y;
X.innr=@(x,y)x'*y;

clear Z;
Z.copy=@(x)x;
Z.scal=@(alpha,x)alpha*x;
Z.zero=@(x)zeros(size(x));
Z.axpy=@(alpha,x,y)alpha*x+y;
Z.innr=@(x,y)x'*y;
Z.prod=@(x,y)x.*y;
Z.id=@(x)ones(size(x));
Z.linv=@(x,y)y./x;
Z.barr=@(x)sum(log(x));
Z.srch=@(x,y)1./max(1./(-y./x));

% Create a bundle of vector spaces
clear VS;
VS.X = X;
VS.Z = Z;

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

% Simple bound constraints on the variables
%
% h(x,y) = [ x ]
%          [ y ]
%
clear h
h.eval=@(x)x;
h.eval_p=@(x,dx)dx;
h.eval_ps=@(x,dy)dy;
h.eval_pps=@(x,dx,dy)zeros(size(x));

% Create the bundle of functions
clear fns;
fns.f=f;
fns.h=h;

% Set a starting guess
pts.x=[2;2.1];
pts.z=[1;1];

% Test the mex file
x=peopt(VS,fns,pts,'run.peopt');
