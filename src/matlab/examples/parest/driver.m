% This demonstrates how to solve a simple parameter estimation problem
%
% min .5 || x2 - d ||^2 + .5 beta || x1 || ^2 s.t. (sum_i A_i x1_i)x_2 = b
%
% where A_i, d, and b are randomly generated.  Note, since they're randomly
% generated, it's possible to create a degenerate problem that causes
% problems for our linear solvers.

function run()

% Make sure that Optizelle is in the path
addpath('../../');

% Create the vector spaces
clear X1;
X1.copy=@(x)x;
X1.scal=@(alpha,x)alpha*x;
X1.zero=@(x)zeros(size(x));
X1.axpy=@(alpha,x,y)alpha*x+y;
X1.innr=@(x,y)x'*y;

clear X2;
X2.copy=@(x)x;
X2.scal=@(alpha,x)alpha*x;
X2.zero=@(x)zeros(size(x));
X2.axpy=@(alpha,x,y)alpha*x+y;
X2.innr=@(x,y)x'*y;

% Set the size of the problem
m = 5; % Size of x1
n = 7; % Size of x2

% Create the data
d=randn(n,1);

% Create the rhs
b=randn(n,1);

% Create the model
clear A;
for i=1:m
    A{i}=randn(n);
end

% ff is the expanded objecitve function with two pieces x{1} and x{2}
beta = 1e-2; % Regularization parameter
ff.eval=@(x) .5*norm(x{2}-d,2)^2 + .5*beta*norm(x{1},2)^2;
ff.grad_1=@(x) beta*x{1};
ff.grad_2=@(x) x{2}-d;
ff.hessvec_11=@(x,dx) beta*dx;
ff.hessvec_12=@(x,dx) zeros(m,1);
ff.hessvec_21=@(x,dx) zeros(n,1);
ff.hessvec_22=@(x,dx) dx;

% gg is the expanded constraint where we split up the variable into two
% pieces x{1} and x{2}.
gg.eval=@(x) g_eval(A,b,x); 
gg.p_1=@(x,dx1) gp_x1(A,x,dx1); 
gg.p_2=@(x,dx2) gp_x2(A,x,dx2); % Need inverse
gg.p_2_inv=@(x,dx2) gp_x2_inv(A,x,dx2); 
gg.ps_1=@(x,dy) gps_x1(A,x,dy); 
gg.ps_2=@(x,dy) gps_x2(A,x,dy); % Need inverse
gg.ps_2_inv=@(x,dy) gps_x2_inv(A,x,dy); 
gg.pps_11=@(x,dx,dy) zeros(m,1);
gg.pps_21=@(x,dx,dy) gps_x2(A,{dx,zeros(n,1)},dy); 
gg.pps_12=@(x,dx,dy) gps_x1(A,{zeros(m,1),dx},dy); 
gg.pps_22=@(x,dx,dy) zeros(n,1);

% Create the equality constrained problem
[X1xX2 f g]=genEqualityConstrained(X1,X2,ff,gg);

% Create a bundle of vector spaces
clear VS;
VS.X = X1xX2;
VS.Y = X2;

% Create the bundle of functions
clear fns;
fns.f=f;
fns.g=g;

% Set a starting guess
pts.x={randn(m,1),randn(n,1)};
pts.dx={randn(m,1),randn(n,1)};
pts.dxx={randn(m,1),randn(n,1)};
pts.y=randn(n,1);
pts.dy=randn(n,1);

% Test the full-space functions 
fprintf('---------Diagnostics on the equality constrained problem---------\n');
Optizelle(VS,fns,pts);

% Solve the full-space problem 
fprintf('\n------------Solving the full-space problem------------\n');
sol=Optizelle(VS,fns,pts,'parest.param');
fprintf('\nWe converged due to: %s\n',sol.opt_stop);

% Form the solution operator.  The user must specify this
phi.eval=@(x)phi_eval(A,b,x);

% Create the unconstrained problem
f = genUnconstrained(X1,X2,ff,gg,phi);

% Create a bundle of vector spaces
clear VS;
VS.X = X1;

% Create the bundle of functions
clear fns;
fns.f=f;

% Set a starting guess
pts.x=randn(m,1);
pts.dx=randn(m,1);
pts.dxx=randn(m,1);
pts.dy=randn(n,1);

% Test the unconstrained functions 
fprintf('\n------------Diagnostics on the unconstrained problem------------\n');
Optizelle(VS,fns,pts);

% Solve the reduced-space problem 
fprintf('\n------------Solving the reduced-space problem------------\n');
sol=Optizelle(VS,fns,pts,'parest.param');
fprintf('\nWe converged due to: %s\n',sol.opt_stop);
end

% Create the functions for the constraint
function z = g_eval(A,b,x)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x{1}(i)*A{i};
    end

    % Now, find (sum A_i x1_i)x2 - f
    z=B*x{2}-b;
end

function z = gp_x1(A,x,dx1)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i dx1_i
    B = zeros(n);
    for i=1:m
        B = B + dx1(i)*A{i};
    end

    % Now, find (sum A_i dx1_i)y
    z=B*x{2};
end

function z = gp_x2(A,x,dx2)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x{1}(i)*A{i};
    end

    % Now, find (sum A_i x1_i)dx2
    z=B*dx2;
end

function z = gp_x2_inv(A,x,dx2)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x{1}(i)*A{i};
    end

    % Now, solve (sum A_i x1_i) z = dx2
    z=B\dx2;
end

function z = gps_x1(A,x,dy)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);

    % Form the adjoint
    z=zeros(m,1);
    for i=1:m
        z(i)=x{2}'*A{i}'*dy;
    end
end

function z = gps_x2(A,x,dy)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x{1}(i)*A{i}';
    end

    % Now, find (sum A_i' x1_i)dy
    z=B*dy;
end

function z = gps_x2_inv(A,x,dy)
    % Get the sizes
    m = size(x{1},1);
    n = size(x{2},1);
    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x{1}(i)*A{i}';
    end

    % Now, solve (sum A_i' x1_i) z = dy
    z=B\dy;
end

% Create the functions for the solution operator
function z = phi_eval(A,b,x)
    % Cache the result if possible
    persistent zz;
    persistent xx;

    % Get the sizes
    m = size(x,1);
    n = length(b);

    % Check if we've already evaluated the function.
    if length(xx)==m
        if norm(x-xx)/(1e-16+norm(x)) < 1e-15
            z=zz;
            return;
        end
    end

    
    % First, form sum A_i x1_i
    B = zeros(n);
    for i=1:m
        B = B + x(i)*A{i};
    end

    % Now, solve (sum A_i x_i) z = b 
    z=B\b;

    % Cache the result
    xx=x;
    zz=z;
end
