% Generates the necessary functions for a reduced space formulation of a
% simple, steady-state  convection diffusion problem
function [f PH phi] = generate_reduced_space(params)
    % Objective function
    f.eval = @(x)obj_eval(params,x);
    f.grad = @(x)obj_grad(params,x);
    f.hessvec = @(x,dx)obj_hv(params,x,dx);

    % Preconditioner for the objective function
    PH.eval = @(state,dx)obj_hv_inv(params,state.x,dx);

    % Finds the state solution
    phi = @(x)state_uncached(params,x,rhs(params,x));
end

%---ObjectiveGradient0---
% Evaluates the objective
function z = obj_eval(params,x)
    % Cached objective evaluation.  Really, this only saves us the first
    % objective evaluation as the subsequent evaluations are cached by
    % Optizelle
    global ocache

    % Performance diagnostis
    global diagnostics

    % Grab the cached objective evaluation when possible
    if ~isempty(ocache)  && isequal(x,ocache.x)
        z = ocache.eval;
        diagnostics.used_cached_objective = diagnostics.used_cached_objective+1;
    else
        % We don't use the caching state solve here because the objective
        % may be evaluated at multiple points during a single optimization
        % iteration, primarily for globalization.  This differs from the
        % gradient and Hessian-vector product, which are both evalated at a
        % fixed point each iteration.
        u = state_uncached(params,x,rhs(params,x));

        % Evaluate the objective
        z = 0.5 * norm(u-params.d)^2;
    end
end

% Evaluates the gradient
function grad = obj_grad(params,x)
    % Cached objective evaluation
    global ocache

    % Solve for the current solution
    u = state(params,x,rhs(params,x));

    % Cached the state solution globally for the objective
    if isempty(ocache) || ~isequal(x,ocache.x)
        ocache.x = x;
        ocache.eval = 0.5 * norm(u-params.d)^2;
    end

    % Set each element of the gradient
    grad = zeros(2,1);
    for i=1:2
        grad(i) = innr( ...
            u-params.d, ...
            -state(params,x,op_p(i,params,x)*u - rhs_p(i,params,x)));
    end
end
%---ObjectiveGradient1---

%---Hessian0---
% Evaluates the Hessian-vector product
function hv = obj_hv(params,x,dx)
    hv = hessian(params,x)*dx;
end

% Finds the Hessian
function H = hessian(params,x)
    % Keep track of where the construction occurs
    persistent cache

    % Performance diagnostics
    global diagnostics

    % Cache the Hessian when required
    if isempty(cache) || ~isequal(x,cache.x)
        % Save the point we're evaluating the Hessian at
        cache.x = x;

        % Solve for the current solution
        u = state(params,x,rhs(params,x));

        % Calculate the Hessian
        cache.H = zeros(2);
        innr = @(x,y)x'*y;
        for j=1:2
            for i=1:j
                cache.H(i,j) = ...
                    innr( ...
                        -state(params,x, ...
                            op_p(j,params,x)*u - rhs_p(j,params,x)), ...
                        -state(params,x, ...
                            op_p(i,params,x)*u - rhs_p(i,params,x))) + ...
                    innr(u-params.d, ...
                        state(params,x, ...
                            op_p(j,params,x) * state(params,x, ...
                                op_p(i,params,x)*u-rhs_p(i,params,x))))+ ...
                    innr(u-params.d, ...
                        state(params,x, ...
                            op_p(i,params,x) * state(params,x, ...
                                op_p(j,params,x) * u))) + ...
                    innr(u-params.d, ...
                        -state(params,x, ...
                            op_p(i,params,x) * ...
                                state(params,x,rhs_p(j,params,x))));

            end
        end
        cache.H(2,1)=cache.H(1,2);

        % Keep track that we cache a Hessian
        diagnostics.hessian_cached = diagnostics.hessian_cached+1;
    end

    % Evaluate the Hessian-vector product
    H = cache.H;
end
%---Hessian1---

%---HessianInv0---
% Evaluates the inverse of the Hessian applied to a vector
function ihv = obj_hv_inv(params,x,dx)
    % Keep track of where the factorization occurs
    persistent cache

    % Performance diagnostics
    global diagnostics

    % Cache the Hessian factorization when required
    if isempty(cache) || ~isequal(x,cache.x)
        % Save the point we're factorizing the Hessian factorization at
        cache.x = x;

        % Grab the current Hessian
        H = hessian(params,x);

        % Factorize the Hessian
        [cache.l cache.u cache.p]=lu(H,'vector');

        % Keep track that we cache a Hessian factorization
        diagnostics.hessian_factorization_cached = ...
            diagnostics.hessian_factorization_cached+1;
    end

    % Apply the inverse to the direction
    ihv = cache.u\(cache.l\dx(cache.p));
end
%---HessianInv1---

% Adjust the forcing function with the pieces for the boundary conditions
function z = rhs(params,x)
    z = params.f - x(1)*params.Ahat - x(2)*params.Bhat;
end

% Finds the derivative of the rhs function adjusted with the boundary
% conditions
function z = rhs_p(which,params,x)
    if which==1
        z = -params.Ahat;
    else
        z = -params.Bhat;
    end
end

% Grabs the differential operator
function z = op(params,x)
    z = x(1)*params.A+x(2)*params.B;
end

% Grabs the derivative of the differential operator
function z = op_p(which,params,x)
    if which==1
        z = params.A;
    else
        z = params.B;
    end
end

%---StateSolve0---
% Solves the discretized PDE with caching
function z = state(params,x,rhs)
    % Keep track of where the solve occurs
    persistent cache

    % Performance diagnostics
    global diagnostics

    % Cache the factorization when required
    if isempty(cache) || ~isequal(x,cache.x)
        % Save the point we're factorizing at
        cache.x = x;

        % Factorize the operator
        [cache.l cache.u cache.p cache.q cache.r] = ...
            lu(op(params,x),'vector');

        % Keep track that we did a new factorization
        diagnostics.state_factorization_cached = ...
            diagnostics.state_factorization_cached+1;
    end

    % Solve the linear system
    z = zeros(size(rhs));
    z(cache.q) = cache.u\(cache.l\(cache.r(:,cache.p)\rhs));
end
%---StateSolve1---

% Solves the discretized PDE without caching
function z = state_uncached(params,x,rhs)
    z = op(params,x)\rhs;
end

% Inner product for Rm
function z = innr(x,y)
    z=x'*y;
end
