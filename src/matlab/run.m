% This tests our ability to create a use a vector space through a mex file

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
Z.barr=@(x,y)sum(log(x));
Z.srch=@(x,y)1./max(1./(-x./y));

% Create a bundle of vector spaces
clear VS;
VS.X = X;
VS.Y = Y;
VS.Z = Z;

% Create the functions
clear f;
f.eval=@(x)(1-x(1))^2*100*(x(2)-x(1)^2)^2;
f.grad=@(x) ...
    [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); ...
    200*(x(2)-x(1)^2)];
f.hessvec=@(x,dx) ...
    [(1200*sq(x(1))-400*x(2)+2)*dx(1)-400*x(1)*dx(2); ...
    -400*x(1)*dx(1) + 200*dx(2)];

clear g;
g.eval=@(x) ...
    [x(1)+2.*x(2)-1.; ...
    2.*x(1)+x(2)-1.];




% Create the bundle of functions
fns.f=f;

% Test the mex file
getProblemClass(VS,fns);
