% This demonstrates how to run a finite difference check on both the
% objective functions as well as the constraints

% Make sure that peopt is in the path
addpath('../../');

% Create the vector spaces
clear X;
X.init=@(x)x;
X.copy=@(x)x;
X.scal=@(alpha,x)alpha*x;
X.zero=@(x)zeros(size(x));
X.axpy=@(alpha,x,y)alpha*x+y;
X.innr=@(x,y)x'*y;

clear Y;
Y.init=@(x)x;
Y.copy=@(x)x;
Y.scal=@(alpha,x)alpha*x;
Y.zero=@(x)zeros(size(x));
Y.axpy=@(alpha,x,y)alpha*x+y;
Y.innr=@(x,y)x'*y;

clear Z;
Z.init=@(x)x;
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
VS.Y = Y;
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


% Define some utility function where
%
% g(x,y)= [ cos(x)sin(y) ]
%         [ 3 x^2 y + y^3]
%         [ log(x) + 3y^5]
%
clear g;
g.eval=@(x) ...
    [cos(x(1))*sin(x(2)); ...
    3.*(x(1)^2)*x(2)+(x(2)^3); ...
    log(x(1))+3.*(x(2))^5];
g.eval_p=@(x,dx) ...
    [-sin(x(1))*sin(x(2))*dx(1)+ ...
        cos(x(1))*cos(x(2))*dx(2); ...
    6.*x(1)*x(2)*dx(1)+ ...
        (3.*(x(1)^2)+3.*(x(2)^2))*dx(2); ...
    1./x(1)*dx(1)+ ...
        15.*(x(2)^4)*dx(2)];
g.eval_ps=@(x,dy) ...
    [-sin(x(1))*sin(x(2))*dy(1)+ ...
          6.*x(1)*x(2)*dy(2)+ ...
          1./x(1)*dy(3); ...
    cos(x(1))*cos(x(2))*dy(1)+ ...
          (3.*(x(1)^2)+3.*(x(2)^2))*dy(2)+ ...
          15.*(x(2)^4)*dy(3)];
g.eval_pps=@(x,dx,dy) ...
    [(-cos(x(1))*dx(1)*sin(x(2))-sin(x(1))*cos(x(2))*dx(2))*dy(1)+ ...
        (6.*dx(1)*x(2) + 6.*x(1)*dx(2))*dy(2)+ ...
        (-1./(x(1)^2)*dx(1))*dy(3);
    (-sin(x(1))*dx(1)*cos(x(2))-cos(x(1))*sin(x(2))*dx(2))*dy(1)+ ...
        (6.*x(1)*dx(1)+6.*x(2)*dx(2))*dy(2)+ ...
        (60.*(x(2)^3)*dx(2))*dy(3)];

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
fns.g=g;
fns.h=h;

% Set a starting guess
pts.x=[2;2.1];
pts.dx=randn(2,1);
pts.dxx=randn(2,1);
pts.dy=randn(3,1);
pts.dz=randn(2,1);

% Test the mex file
peopt(VS,fns,pts);
