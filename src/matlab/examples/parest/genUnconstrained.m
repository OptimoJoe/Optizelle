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
% f : The reduced objective function .
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

% This function needs to cache both x and dx since it's called in more than one
% place.
% g'_2(x,phi(x)) phi'(x) = -g'_1(x,phi(x)) dx
phi.p=@(x,dx)p(X1,x,dx,gg,phi,kry_tol,kry_rst,kry_iter_max);
% phi'(x)*dy = -g'_1(x,phi(x))* inv(g'_2(x,phi(x))*) dy
phi.ps=@(x,dy) ...
    feval(@(xx) ...
        X1.scal(-1., ...
            gg.ps_1(xx, ...
                mygmres( ...
                    @(dy)gg.ps_2(xx,dy), ...
                    dy, ... 
                    kry_rst,kry_tol,kry_iter_max, ...
                    @(dy)gg.ps_2_inv(xx,dy)))), ...
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
phi.pps=@(x,dx,dy_hat) ...
    X1.scal(-1., ...
        feval(@(xx) ...
            feval(@(dy) ...
                feval(@(phi_p_x_dx) ...
                    X1.axpy(1.,gg.pps_11(xx,dx,dy), ...
                    X1.axpy(1.,gg.pps_12(xx,phi_p_x_dx,dy),...
                    phi.ps(x, ...
                        X2.axpy(1., ...
                            gg.pps_21(xx,dx,dy), ...
                            gg.pps_22(xx,phi_p_x_dx,dy))))), ...
                    phi.p(x,dx)), ...
                mygmres( ...
                    @(dy)gg.ps_2(xx,dy), ...
                    dy_hat, ...
                    kry_rst,kry_tol,kry_iter_max, ...
                    @(dy)gg.ps_2_inv(xx,dy))), ...
            {x,phi.eval(x)}));

% Form the objective function
f.eval=@(x)ff.eval({x,phi.eval(x)});
f.grad=@(x) ...
    feval(@(xx) ...
        X1.axpy(1., ...
            ff.grad_1(xx), ...
            phi.ps(x,ff.grad_2(xx))), ...
        {x,phi.eval(x)});
f.hessvec=@(x,dx) ...
    feval(@(xx) ...
        feval(@(phi_p_x_dx) ...
            X1.axpy(1.,ff.hessvec_11(xx,dx), ...
            X1.axpy(1.,ff.hessvec_12(xx,phi_p_x_dx), ...
            X1.axpy(1.,phi.pps(x,dx,ff.grad_2(xx)), ...
                       phi.ps(x, ...
                            X2.axpy(1., ...
                                ff.hessvec_21(xx,dx),...
                                ff.hessvec_22(xx,phi_p_x_dx)))))), ...
            phi.p(x,dx)), ...
        {x,phi.eval(x)});

end

% Derivative of the solution operator 
function z=p(X1,x,dx,gg,phi,kry_tol,kry_rst,kry_iter_max)
    % Cache the result if possible
    persistent zz;
    persistent xx;
    persistent dxx;

    % Check if we've already evaluated the function.
    if ~isempty(xx)
        res_x = X1.axpy(-1.,xx,x);
        res_dx = X1.axpy(-1.,dxx,dx);
        if( sqrt(X1.innr(res_x,res_x))/(1e-16+sqrt(X1.innr(x,x))) < 1e-15 && ...
            sqrt(X1.innr(res_dx,res_dx))/(1e-16+sqrt(X1.innr(dx,dx))) < 1e-15)
            z=zz;
            return;
        end
    end

    % Otherwise, calculate the result
    % g'_2(x,phi(x)) phi'(x) = -g'_1(x,phi(x)) dx
    z=feval(@(xx) ...
        mygmres( ...
            @(dx)gg.p_2(xx,dx), ...
            X1.scal(-1.,gg.p_1(xx,dx)), ...
            kry_rst,kry_tol,kry_iter_max, ...
            @(dx)gg.p_2_inv(xx,dx)), ...
        {x,phi.eval(x)});

    % Now, cache everything
    xx=x;
    dxx=dx;
    zz=z;
end

% A function to eliminate the messages from MATLAB's gmres routines
function z = mygmres(A,b,kry_rst,kry_tol,kry_iter_max,M)
    [z flag]=gmres(A,b,kry_rst,kry_tol,kry_iter_max,M);
end
