% Given an equality constrained problem defined on two variables,
%
% min_{x1,x2} f(x1,x2) st g(x1,x2)=0,
%
% we can define the derivatives of these functions in several parts, such as
% grad_1 f(x1,x2), grad_2 f(x1,x2), and so on.  From these pieces, this function
% produces the composite functions for f and g, which are then appropriate for
% equality constrained optimization.
%
% function [X1xX2 f g]=genEqualityConstrained(X1,X2,ff,gg)
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
%
% OUTPUTS
%
% X1xX2: The composite vector space and operators for X1 and X2.  The variables
%        in the composite vector space are simply cell arrays where x{1} and
%        x{2} correspond to the returned variables.
% f : The composite objective function.
% g : The composite equality constraints.
%

function [X1xX2 f g]=genEqualityConstrained(X1,X2,ff,gg)

% Generate the composite vector space
X1xX2.copy=@(x){X1.copy(x{1}),X2.copy(x{2})};
X1xX2.scal=@(alpha,x){X1.scal(alpha,x{1}),X2.scal(alpha,x{2})};
X1xX2.zero=@(x){X1.zero(x{1}),X2.zero(x{2})};
X1xX2.axpy=@(alpha,x,y){X1.axpy(alpha,x{1},y{1}),X2.axpy(alpha,x{2},y{2})};
X1xX2.innr=@(x,y)X1.innr(x{1},y{1})+X2.innr(x{2},y{2});

% f is the composite objective function built from ff.
f.eval=@(x)ff.eval(x);
f.grad=@(x){ff.grad_1(x),ff.grad_2(x)};
f.hessvec=@(x,dx){ ...
    ff.hessvec_11(x,dx{1})+ff.hessvec_12(x,dx{2}), ...
    ff.hessvec_21(x,dx{1})+ff.hessvec_22(x,dx{2})};

% g is the composite constraint built from gg. 
g.eval=@(x)gg.eval(x);
g.p=@(x,dx)gg.p_1(x,dx{1}) + gg.p_2(x,dx{2});
g.ps=@(x,dy){gg.ps_1(x,dy),gg.ps_2(x,dy)};
g.pps=@(x,dx,dy)X1xX2.axpy(1.0, ...
            {gg.pps_11(x,dx{1},dy),gg.pps_21(x,dx{1},dy)}, ...
            {gg.pps_12(x,dx{2},dy),gg.pps_22(x,dx{2},dy)});
