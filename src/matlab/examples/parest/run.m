% This demonstrates how to do simple parameter estimation.
function run()

% Make sure that peopt is in the path
addpath('../../');

% Create the vector spaces
clear X1;
X1.init=@(x)x;
X1.copy=@(x)x;
X1.scal=@(alpha,x)alpha*x;
X1.zero=@(x)zeros(size(x));
X1.axpy=@(alpha,x,y)alpha*x+y;
X1.innr=@(x,y)x'*y;

clear X2;
X2.init=@(x)x;
X2.copy=@(x)x;
X2.scal=@(alpha,x)alpha*x;
X2.zero=@(x)zeros(size(x));
X2.axpy=@(alpha,x,y)alpha*x+y;
X2.innr=@(x,y)x'*y;

clear X1xX2;
X1xX2.init=@(x){X1.init(x{1}),X2.init(x{2})};
X1xX2.copy=@(x){X1.copy(x{1}),X2.copy(x{2})};
X1xX2.scal=@(alpha,x){X1.scal(alpha,x{1}),X2.scal(alpha,x{2})};
X1xX2.zero=@(x){X1.zero(x{1}),X2.zero(x{2})};
X1xX2.axpy=@(alpha,x,y){X1.axpy(alpha,x{1},y{1}),X2.axpy(alpha,x{2},y{2})};
X1xX2.innr=@(x,y)X1.innr(x{1},y{1})+X2.innr(x{2},y{2});

% Set the size of the problem

% Size of x1
m = 5;

% Size of x2
n = 7;

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

% f is the compresssed objective function built from the above pieces
f.eval=@(x)ff.eval(x);
f.grad=@(x){ff.grad_1(x),ff.grad_2(x)};
f.hessvec=@(x,dx){ ...
    ff.hessvec_11(x,dx{1})+ff.hessvec_12(x,dx{2}), ...
    ff.hessvec_21(x,dx{1})+ff.hessvec_22(x,dx{2})};

% gg is the expanded constraint where we split up the variable into two
% pieces x{1} and x{2}.
gg.eval=@(x) g_eval(A,b,x); 
gg.eval_p_1=@(x,dx1) gp_x1(A,x,dx1); 
gg.eval_p_2=@(x,dx2) gp_x2(A,x,dx2); % Need inverse
gg.eval_p_2_inv=@(x,dx2) gp_x2_inv(A,x,dx2); 
gg.eval_ps_1=@(x,dy) gps_x1(A,x,dy); 
gg.eval_ps_2=@(x,dy) gps_x2(A,x,dy); % Need inverse
gg.eval_ps_2_inv=@(x,dy) gps_x2_inv(A,x,dy); 
gg.eval_pps_11=@(x,dx,dy) zeros(m,1);
gg.eval_pps_21=@(x,dx,dy) gps_x2(A,{dx,zeros(n,1)},dy); 
gg.eval_pps_12=@(x,dx,dy) gps_x1(A,{zeros(m,1),dx},dy); 
gg.eval_pps_22=@(x,dx,dy) zeros(n,1);

% g is the compressed constraint that peopt sees.  This combines the 
% expanded functions from above.
g.eval=@(x)gg.eval(x);
g.eval_p=@(x,dx)gg.eval_p_1(x,dx{1}) + gg.eval_p_2(x,dx{2});
g.eval_ps=@(x,dy){gg.eval_ps_1(x,dy),gg.eval_ps_2(x,dy)};
g.eval_pps=@(x,dx,dy)X1xX2.axpy(1.0, ...
            {gg.eval_pps_11(x,dx{1},dy),gg.eval_pps_21(x,dx{1},dy)}, ...
            {gg.eval_pps_12(x,dx{2},dy),gg.eval_pps_22(x,dx{2},dy)});

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
pts.dy=randn(n,1);

% Test the full-space functions 
fprintf('--------------Diagnostics on the full-space problem--------------\n');
peopt(VS,fns,pts);


% Form the solution operator
kry_tol = 1e-10;
kry_rst = n;
kry_iter_max = n;
% The user must specify this
phi.eval=@(x)phi_eval(A,b,x);
% g'_2(x,phi(x)) phi'(x) = -g'_1(x,phi(x)) dx
phi.eval_p=@(x,dx) ...
    feval(@(xx) ...
        mygmres( ...
            @(dx)gg.eval_p_2(xx,dx), ...
            X1.scal(-1.,gg.eval_p_1(xx,dx)), ...
            kry_rst,kry_tol,kry_iter_max, ...
            @(dx)gg.eval_p_2_inv(xx,dx)), ...
        {x,phi.eval(x)});
% phi'(x)*dy = -g'_1(x,phi(x))* inv(g'_2(x,phi(x))*) dy
phi.eval_ps=@(x,dy) ...
    feval(@(xx) ...
        X1.scal(-1., ...
            gg.eval_ps_1(xx, ...
                mygmres( ...
                    @(dy)gg.eval_ps_2(xx,dy), ...
                    dy, ... 
                    kry_rst,kry_tol,kry_iter_max, ...
                    @(dy)gg.eval_ps_2_inv(xx,dy)))), ...
        {x,phi.eval(x)});
% (phi''(x)dx)*dy_hat = 
% -(   (g''_11(x,phi(x))dx)*dy
%    + (g''_12(x,phi(x))phi'(x)dx)*dy
%    + phi'(x)* (g''_21(x,phi(x))dx)*dy
%                + g''_22(x,phi(x))phi'(x)dx)*dy)
%  ) dy
%
% where we solve for dy in g'_2(x,phi(x))* dy = dy_hat.
%
phi.eval_pps=@(x,dx,dy_hat) ...
    X1.scal(-1., ...
        feval(@(xx) ...
            feval(@(dy) ...
                feval(@(phi_p_x_dx) ...
                    X1.axpy(1.,gg.eval_pps_11(xx,dx,dy), ...
                    X1.axpy(1.,gg.eval_pps_12(xx,phi_p_x_dx,dy),...
                    phi.eval_ps(x, ...
                        X2.axpy(1., ...
                            gg.eval_pps_21(xx,dx,dy), ...
                            gg.eval_pps_22(xx,phi_p_x_dx,dy))))), ...
                    phi.eval_p(x,dx)), ...
                mygmres( ...
                    @(dy)gg.eval_ps_2(xx,dy), ...
                    dy_hat, ...
                    kry_rst,kry_tol,kry_iter_max, ...
                    @(dy)gg.eval_ps_2_inv(xx,dy))), ...
            {x,phi.eval(x)}));

% Form the objective function
theta.eval=@(x)ff.eval({x,phi.eval(x)});
theta.grad=@(x) ...
    feval(@(xx) ...
        X1.axpy(1., ...
            ff.grad_1(xx), ...
            phi.eval_ps(x,ff.grad_2(xx))), ...
        {x,phi.eval(x)});
theta.hessvec=@(x,dx) ...
    feval(@(xx) ...
        feval(@(phi_p_x_dx) ...
            X1.axpy(1.,ff.hessvec_11(xx,dx), ...
            X1.axpy(1.,ff.hessvec_12(xx,phi_p_x_dx), ...
            X1.axpy(1.,phi.eval_pps(x,dx,ff.grad_2(xx)), ...
                       phi.eval_ps(x, ...
                            X2.axpy(1., ...
                                ff.hessvec_21(xx,dx),...
                                ff.hessvec_22(xx,phi_p_x_dx)))))), ...
            phi.eval_p(x,dx)), ...
        {x,phi.eval(x)});

% Create a bundle of vector spaces
clear VS;
VS.X = X1;
VS.Y = X2;

% Create the bundle of functions
clear fns;
fns.f=theta;
fns.g=phi;

% Set a starting guess
pts.x=randn(m,1);
pts.dx=randn(m,1);
pts.dxx=randn(m,1);
pts.dy=randn(n,1);

% Test the reduced-space functions 
fprintf('\n------------Diagnostics on the reduced-space problem------------\n');
peopt(VS,fns,pts);

% Create a bundle of vector spaces
clear VS;
VS.X = X1;

% Create the bundle of functions
clear fns;
fns.f=theta;

% Set a starting guess
pts.x=randn(m,1);
pts.dx=randn(m,1);
pts.dxx=randn(m,1);

% Solve the reduced-space problem 
fprintf('\n------------Solving the reduced-space problem------------\n');
%x=peopt(VS,fns,pts,'run.peopt');
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

% A function to eliminate the messages from MATLAB's gmres routines
function z = mygmres(A,b,kry_rst,kry_tol,kry_iter_max,M)
    [z flag]=gmres(A,b,kry_rst,kry_tol,kry_iter_max,M);
end
