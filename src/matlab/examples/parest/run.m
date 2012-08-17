% This demonstrates how to do simple parameter estimation.
function run()

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

clear XxY;
XxY.init=@(x)x;
XxY.copy=@(x)x;
XxY.scal=@(alpha,x){alpha*x{1},alpha*x{2}};
XxY.zero=@(x){zeros(size(x{1})),zeros(size(x{2}))};
XxY.axpy=@(alpha,x,y){alpha*x{1}+y{1},alpha*x{2}+y{2}};
XxY.innr=@(x,y)x{1}'*y{1}+x{2}'*y{2};

% Set the size of the problem

% Size of x1
m = 5;

% Size of x2
n = 7;

% Create the data
d=randn(n,1);

% Create the model
clear A;
for i=1:m
    A{i}=randn(n);
end

% ff is the expanded objecitve function with two pieces x{1} and x{2}
ff.eval=@(x) .5*norm(x{2}-d,2)^2;
ff.grad_1=@(x) zeros(m,1); 
ff.grad_2=@(x) x{2}-d;
ff.hessvec_11=@(x,dx) zeros(m,1);
ff.hessvec_12=@(x,dx) zeros(m,1);
ff.hessvec_21=@(x,dx) zeros(n,1);
ff.hessvec_22=@(x,dx) dx;

% f is the compresssed objective function built from the above pieces
f.eval=@(x)ff.eval(x);
f.grad=@(x){ff.grad_1(x),ff.grad_2(x)};
f.hessvec=@(x,dx){ ...
    ff.hessvec_11(x,dx{1})+ff.hessvec_12(x,dx{2}), ...
    ff.hessvec_21(x,dx{1})+ff.hessvec_22(x,dx{2})};

% gg is the expanded constraint where we split up the variable into two
% pieces x{1} and x{2}.
gg.eval=@(x) g_eval(A,x{1},x{2}); 
gg.eval_p_1=@(x,dx1) gp_x(A,x{1},x{2},dx1); 
gg.eval_p_2=@(x,dx2) gp_y(A,x{1},x{2},dx2);
gg.eval_ps_1=@(x,dy) gps_x(A,x{1},x{2},dy); 
gg.eval_ps_2=@(x,dy) gps_y(A,x{1},x{2},dy);
gg.eval_pps_11=@(x,dx,dy)  {zeros(m,1),zeros(n,1)};
gg.eval_pps_12=@(x,dx,dy) {zeros(m,1),zeros(n,1)}; 
gg.eval_pps_21=@(x,dx,dy) {zeros(m,1),zeros(n,1)};
gg.eval_pps_22=@(x,dx,dy) {zeros(m,1),zeros(n,1)};

% g is the compressed constraint that peopt sees.  This combines the 
% expanded functions from above.
g.eval=@(x)gg.eval(x)
g.eval_p=@(x,dx)gg.eval_p_1(x,dx{1}) + gg.eval_p_2(x,dx{2});
g.eval_ps=@(x,dy){gg.eval_ps_1(x,dy),gg.eval_ps_2(x,dy)};
g.eval_pps=@(x,dx,dy)XxY.axpy(1.0, ...
            gg.eval_pps_11(x,dx{1},dy), ...
            XxY.axpy(1.0, ...
                 gg.eval_pps_12(x,dx{2},dy), ...
                 XxY.axpy(1.0, ...
                     gg.eval_pps_21(x,dx{1},dy), ...
                     gg.eval_pps_22(x,dx{2},dy))));


% Create a bundle of vector spaces
clear VS;
VS.X = XxY;
VS.Y = Y;

% Create the bundle of functions
clear fns;
fns.f=f;
fns.g=g;

% Set a starting guess
pts.x={randn(m,1),randn(n,1)};
pts.dx={randn(m,1),randn(n,1)};
pts.dxx={randn(m,1),randn(n,1)};
pts.dy=randn(n,1);

% Test the mex file
peopt(VS,fns,pts);
end

% Create the functions for the constraint
function z = g_eval(A,x,y)
    % Get the sizes
    m = size(x,1);
    n = size(y,1);
    
    % First, form sum A_i x_i
    B = x(1)*A{1};
    for i=2:m
        B = B + x(i)*A{i};
    end

    % Now, find (sum A_i x_i)y
    z=B*y;
end

function z = gp_x(A,x,y,dx)
    % Get the sizes
    m = size(x,1);
    n = size(y,1);
    
    % First, form sum A_i dx_i
    B = x(1)*A{1};
    for i=2:m
        B = B + dx(i)*A{i};
    end

    % Now, find (sum A_i dx_i)y
    z=B*y;
end

function z = gp_y(A,x,y,dy)
    % Get the sizes
    m = size(x,1);
    n = size(y,1);
    
    % First, form sum A_i x_i
    B = x(1)*A{1};
    for i=2:m
        B = B + x(i)*A{i};
    end

    % Now, find (sum A_i x_i)dy
    z=B*dy;
end

function z = gps_x(A,x,y,dz)
    % Get the sizes
    m = size(x,1);
    n = size(y,1);

    % Form the adjoint
    z=zeros(m,1);
    for i=1:m
        z(i)=y'*A{i}*dz;
    end
end

function z = gps_y(A,x,y,dz)
    % Get the sizes
    m = size(x,1);
    n = size(y,1);
    
    % First, form sum A_i' x_i
    B = x(1)*A{1}';
    for i=2:m
        B = B + x(i)*A{i}';
    end

    % Now, find (sum A_i' x_i)dy
    z=B*dz;
end
