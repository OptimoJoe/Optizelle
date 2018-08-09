% Given an equality constrained problem defined on two variables,
%
% min_{x1,x2} f(x1,x2) st g(x1,x2)=0,
%
% we can define the derivatives of these functions in several parts, such as
% grad_1 f(x1,x2), grad_2 f(x1,x2), and so on.  Next, given a solution operator
% phi(x1) such that g(x1,phi(x1))=0, we can generate a reduced objective
% function f(x1,phi(x1)), which is then appropriate for unconstrained
% optimization methods.  This is known as a null-space method or a reduced
% space formulation.
%
% function f=genUnconstrained(X1,X2,ff,gg,phi,kry_tol,kry_rst,kry_iter_max)
%
% INPUTS
%
% X1: Vector space for the variable x1.  This is a struct array that contains
%     the elements copy, scal, zero, axpy, and innr
% X2: Vector space for the variable x2.
% ff: Definition of the objective.  This is a struct array that contains
%     the elements eval, grad_1, grad_2, hessvec_11, hessvec_12, hessvec_21,
%     and hessvec_22.
% gg: Definition of the constraint.  This is a struct array that contains
%     the elements eval, p_1, p_2, p_2_inv, ps_1,
%     ps_2, ps_2_inv, pps_11, pps_21, pps_12, and
%     pps_22.  Note, for the inverse operators, only approximate inverses
%     are required.
% phi: Definition of the solution operator.  This is a struct array that
%     constains the element eval.  Note, it is *very* important that
%     phi caches its input.  In other words, the following code will call
%     phi(x1) repeatedly where x1 does not change.  Since the computation
%     of phi can be expensive, it is best that the function phi caches the
%     result and returns this when possible.  Due to the nature of the
%     algorithms involved, we don't generally need to cache more than a
%     single x1 and phi(x1).  Namely, we just need to cache the last pair.
% kry_tol (optional) : GMRES stopping tolernace for the iterative solves
%     required for the algorithm.  Note, we need relatively accurate solves
%     for this method to work.  By default, we use 1e-12.
% kry_rst (optional) : How often GMRES restarts.  By default, we use MATLAB's
%     defaults.
% kry_iter_max (optional) : Maximum number of GMRES iterations.  By default,
%     we use MATLAB's defaults.
%
% OUTPUTS
%
% f : The reduced objective function
%

function f=genUnconstrained(X1,X2,ff,gg,phi,kry_tol,kry_rst,kry_iter_max)

% Figure out our options
if(nargin<5)
    error('Too few arguments to the function.  Check the help.');
elseif(nargin==5)
    kry_tol=1e-12;
    kry_rst=[];
    kry_iter_max=[];
elseif(nargin==6)
    kry_rst=[];
    kry_iter_max=[];
elseif(nargin==7)
    kry_iter_max=[];
elseif(nargin>=9)
    error('Too many arguments to the function.  Check the help.');
end

% Form the solution operator

% phi'(x) dx = - inv(g'_2(x,phi(x))) g'_1(x,phi(x)) dx
phi.p=@(x,dx)p(X1,x,dx,gg,phi,kry_tol,kry_rst,kry_iter_max);

% phi'(x)*dy = -g'_1(x,phi(x))* inv(g'_2(x,phi(x))*) dy
phi.ps=@(x,dy)ps(X1,x,dy,gg,phi,kry_tol,kry_rst,kry_iter_max);

% (phi''(x)dx)*dy =
% -(   (g''_11(x,phi(x))dx)*
%    + (g''_12(x,phi(x))phi'(x)dx)*
%    + phi'(x)* (  (g''_21(x,phi(x))dx)*
%                 +(g''_22(x,phi(x))phi'(x)dx)*)
% ) inv(g'_2(x,phi(x))) dy
phi.pps=@(x,dx,dy) ...
    pps(X1,X2,x,dx,dy,gg,phi,kry_tol,kry_rst,kry_iter_max);

% Form the adjoint solve operator

% psi(x) = inv(g'_2(x,phi(x))*) grad_2 f(x,phi(x))
psi.eval = @(x)psi_eval(X1,X2,x,ff,gg,phi,kry_tol,kry_rst,kry_iter_max);

% Rework the second derivative of the solution operator based on this

% (phi_truncated''(x)dx)* grad_2 f(x,phi(x)) =
% -(   (g''_11(x,phi(x))dx)*
%    + (g''_12(x,phi(x))phi'(x)dx)*
% ) psi(x)
phi.pps_truncated=@(x,dx) ...
    pps_truncated(X1,X2,x,dx,gg,phi,psi,kry_tol,kry_rst,kry_iter_max);

% Form the objective function

% J(x) = f(x,phi(x))
f.eval=@(x)ff.eval({x,phi.eval(x)});

% grad J(x) = grad_1 f(x,phi(x)) + phi'(x)* grad_2 f(x,phi(x))
f.grad_old = @(x)grad_old(X1,x,ff,phi);

% grad J(x) = grad_1 f(x,phi(x)) - g'_1(x,phi(x))* psi(x)
f.grad = @(x)grad(X1,x,ff,gg,phi,psi);

% hess J(x) dx =
%     hess_11 f(x,phi(x)) +
%     hess_12 f(x,phi(x)) phi'(x) +
%     (phi''(x) dx)* grad_2 f(x,phi(x)) +
%     phi'(x)* ( hess_12 f(x,phi(x)) + hess_22 f(x,phi(x)) phi'(x)dx)
f.hessvec_old=@(x,dx) hessvec_old(X1,X2,x,dx,ff,phi);

% hess J(x) dx =
%     hess_11 f(x,phi(x)) +
%     hess_12 f(x,phi(x)) phi'(x) +
%     (phi_truncated''(x) dx)* grad_2 f(x,phi(x)) +
%     phi'(x)* (
%         hess_12 f(x,phi(x)) +
%         hess_22 f(x,phi(x)) phi'(x)dx -
%         (g''_21(x,phi(x))dx)* psi(x) -
%         (g''_22(x,phi(x))phi'(x)dx)* psi(x))
%
% Why go through this mess as opposed to what we had above? If we cache
% properly, we're down to 1 forward and 1 adjoint solve per Krylov iteration.
%
f.hessvec=@(x,dx) hessvec(X1,X2,x,dx,ff,gg,phi,psi);
end

% Derivative of the solution operator.  Requires caching.
function z=p(X1,x,dx,gg,phi,kry_tol,kry_rst,kry_iter_max)
    % Cache the result if possible
    persistent zz;
    persistent xx;
    persistent dxx;

    % Check if we've already evaluated the function.
    if ~isempty(xx) && isequal(xx,x) && isequal(dxx,dx)
        z=zz;
        return;
    end

    % Otherwise, calculate the result

    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    %  z <- phi'(x) dx = - inv(g'_2(x,phi(x))) g'_1(x,phi(x)) dx
    z = mygmres( ...
        @(dx)gg.p_2(x_phix,dx), ...
        X1.scal(-1.,gg.p_1(x_phix,dx)), ...
        kry_rst,kry_tol,kry_iter_max, ...
        @(dx)gg.p_2_inv(x_phix,dx));

    % Now, cache everything
    xx=x;
    dxx=dx;
    zz=z;
end

% Derivative adjoint of the solution operator.  No caching.
function z=ps(X1,x,dy,gg,phi,kry_tol,kry_rst,kry_iter_max)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % z <- phi'(x)*dy = -g'_1(x,phi(x))* inv(g'_2(x,phi(x))*) dy
    z = X1.scal(-1., ...
        gg.ps_1(x_phix, ...
            mygmres( ...
                @(dy)gg.ps_2(x_phix,dy), ...
                dy, ...
                kry_rst,kry_tol,kry_iter_max, ...
                @(dy)gg.ps_2_inv(x_phix,dy))));
end

% Second derivative adjoint of the solution operator.  No caching.
function z=pps(X1,X2,x,dx,dy,gg,phi,kry_tol,kry_rst,kry_iter_max)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % inv_gp2_xphix_dy <- inv(g'_2(x,phi(x))*) dy
    inv_gp2_xphix_dy = mygmres( ...
        @(dy)gg.ps_2(x_phix,dy), ...
        dy, ...
        kry_rst,kry_tol,kry_iter_max, ...
        @(dy)gg.ps_2_inv(x_phix,dy));

    % phi_p_x_dx <- phi'(x)dx
    phi_p_x_dx = phi.p(x,dx);

    % (phi''(x)dx)*dy =
    % -(   (g''_11(x,phi(x))dx)*
    %    + (g''_12(x,phi(x))phi'(x)dx)*
    %    + phi'(x)* (  (g''_21(x,phi(x))dx)*
    %                 +(g''_22(x,phi(x))phi'(x)dx)*)
    % ) inv(g'_2(x,phi(x))) dy
    z = X1.scal(-1., ...
        X1.axpy(1.,gg.pps_11(x_phix,dx,inv_gp2_xphix_dy), ...
        X1.axpy(1.,gg.pps_12(x_phix,phi_p_x_dx,inv_gp2_xphix_dy),...
        phi.ps(x, ...
            X2.axpy(1., ...
                gg.pps_21(x_phix,dx,inv_gp2_xphix_dy), ...
                gg.pps_22(x_phix,phi_p_x_dx,inv_gp2_xphix_dy))))));

end

% Second derivative adjoint of the solution operator with some pieces removed.
% Basically, we're removing them so that we can simplify the Hessian calculation
% later.
function z=pps_truncated(X1,X2,x,dx,gg,phi,psi,kry_tol,kry_rst,kry_iter_max)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % psi_x <- psi(x)
    psi_x = psi.eval(x);

    % phi_p_x_dx <- phi'(x)dx
    phi_p_x_dx = phi.p(x,dx);

    % (phi_truncated''(x)dx)* grad_2 f(x,phi(x)) =
    % -(   (g''_11(x,phi(x))dx)*
    %    + (g''_12(x,phi(x))phi'(x)dx)*
    % ) psi(x)
    z = X1.scal(-1., ...
        X1.axpy(1.,gg.pps_11(x_phix,dx,psi_x), ...
                   gg.pps_12(x_phix,phi_p_x_dx,psi_x)));
end

% Form the adjoint solve operator.  Cached values.
function z=psi_eval(X1,X2,x,ff,gg,phi,kry_tol,kry_rst,kry_iter_max)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % Cache the result if possible
    persistent zz;
    persistent xx;

    % Check if we've already evaluated the function.
    if ~isempty(xx) && isequal(xx,x)
        z=zz;
        return;
    end

    % Otherwise, calculate the result

    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % z <- psi(x) = inv(g'_2(x,phi(x))*) grad_2 f(x,phi(x))
    z = mygmres( ...
        @(dy)gg.ps_2(x_phix,dy), ...
        ff.grad_2(x_phix), ...
        kry_rst,kry_tol,kry_iter_max, ...
        @(dy)gg.ps_2_inv(x_phix,dy));

    % Now, cache everything
    xx=x;
    zz=z;
end

% Gradient of the objective function
function z=grad_old(X1,x,ff,phi)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % z <- grad J(x) = grad_1 f(x,phi(x)) + phi'(x)* grad_2 f(x,phi(x))
    z = X1.axpy(1., ...
        ff.grad_1(x_phix), ...
        phi.ps(x,ff.grad_2(x_phix)));
end

% Gradient of the objective function
function z=grad(X1,x,ff,gg,phi,psi)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % psi_x <- psi(x)
    psi_x = psi.eval(x);

    % z <- grad J(x) = grad_1 f(x,phi(x)) - g'_1(x,phi(x))* psi(x)
    z = X1.axpy(1., ...
        ff.grad_1(x_phix), ...
        X1.scal(-1., ...
            gg.ps_1(x_phix,psi_x)));
end

% Hessian-vector product of the objective function
function z=hessvec_old(X1,X2,x,dx,ff,phi)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % phi_p_x_dx <- phi'(x)dx
    phi_p_x_dx = phi.p(x,dx);

    % z <- hess J(x) dx =
    %     hess_11 f(x,phi(x)) +
    %     hess_12 f(x,phi(x)) phi'(x) +
    %     (phi''(x) dx)* grad_2 f(x,phi(x)) +
    %     phi'(x)* ( hess_12 f(x,phi(x)) + hess_22 f(x,phi(x)) phi'(x)dx)
    z = X1.axpy(1.,ff.hessvec_11(x_phix,dx), ...
        X1.axpy(1.,ff.hessvec_12(x_phix,phi_p_x_dx), ...
        X1.axpy(1.,phi.pps(x,dx,ff.grad_2(x_phix)), ...
            phi.ps(x, ...
                X2.axpy(1., ...
                    ff.hessvec_21(x_phix,dx),...
                    ff.hessvec_22(x_phix,phi_p_x_dx))))));
end

% Hessian-vector product of the objective function
function z=hessvec(X1,X2,x,dx,ff,gg,phi,psi)
    %x_phix <- (x,phi(x))
    x_phix = {x,phi.eval(x)};

    % phi_p_x_dx <- phi'(x)dx
    phi_p_x_dx = phi.p(x,dx);

    % psi_x <- psi(x)
    psi_x = psi.eval(x);

    % z<- hess J(x) dx =
    %     hess_11 f(x,phi(x)) +
    %     hess_12 f(x,phi(x)) phi'(x) +
    %     (phi_truncated''(x) dx)* grad_2 f(x,phi(x)) +
    %     phi'(x)* (
    %         hess_12 f(x,phi(x)) +
    %         hess_22 f(x,phi(x)) phi'(x)dx -
    %         (g''_21(x,phi(x))dx)* psi(x) -
    %         (g''_22(x,phi(x))phi'(x)dx)* psi(x))
    %
    z = X1.axpy(1.,ff.hessvec_11(x_phix,dx), ...
        X1.axpy(1.,ff.hessvec_12(x_phix,phi_p_x_dx), ...
        X1.axpy(1.,phi.pps_truncated(x,dx), ...
                   phi.ps(x, ...
                       X2.axpy(1.,ff.hessvec_21(x_phix,dx),...
                       X2.axpy(1.,ff.hessvec_22(x_phix,phi_p_x_dx), ...
                       X2.scal(-1., ...
                           X2.axpy(1., ...
                               gg.pps_21(x_phix,dx,psi_x), ...
                               gg.pps_22(x_phix,phi_p_x_dx,psi_x)))))))));
end

% A function to eliminate the messages from MATLAB's gmres routines
function z = mygmres(A,b,kry_rst,kry_tol,kry_iter_max,M)
    warning off;
    [z flag]=gmres(A,b,kry_rst,kry_tol,kry_iter_max,M);
    warning on;
end
