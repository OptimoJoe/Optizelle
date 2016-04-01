% Generates the necessary functions for a full space formulation of a simple,
% steady-state  convection diffusion problem
function [f g schur phi] = generate_full_space(params)
    % Objective function
    f.eval = @(x)my_obj(params,x);
    f.grad = @(x)obj_grad(params,x);
    f.hessvec = @(x,dx)obj_hv(params,x,dx);

    % Equality constraint
    g.eval = @(x)eq_eval(params,x);
    g.p = @(x,dx)eq_p(params,x,dx);
    g.ps = @(x,dy)eq_ps(params,x,dy);
    g.pps = @(x,dx,dy)eq_pps(params,x,dx,dy);

    % Preconditioner for the equality constraint
    schur.eval = @(state,dx)eq_schur(params,state.x,dx);
   
    % Finds the state solution 
    phi = @(x)state_uncached(params,x,rhs(params,x));
end

% Evaluates the objective 
function z = my_obj(params,x)
    z = 0.5 * norm(x(params.idx.u)-params.d)^2;
end

% Evaluates the gradient
function grad = obj_grad(params,x)
    grad = zeros(params.nx+2,1);
    grad(params.idx.u) = x(params.idx.u)-params.d;
end

% Evaluates the Hessian-vector product
function hv = obj_hv(params,x,dx)
    hv = zeros(params.nx+2,1);
    hv(params.idx.u) = dx(params.idx.u);
end

% Evaluates the equality constraint
function z = eq_eval(params,x)
    z = op(params,x)*x(params.idx.u)-rhs(params,x);
end

% Evaluates the derivative of the equality constraint 
function z = eq_p(params,x,dx)
    z = deriv(params,x)*dx;
end

% Evaluates the adjoint of the derivative of the equality constraint 
function z = eq_ps(params,x,dy)
    z = deriv(params,x)'*dy;
end

% Evaluates the adjoint of second derivative of the equality constraint
function z = eq_pps(params,x,dx,dy)
    z = deriv2(params,x,dy)*dx;
end

% Evaluates the Schur preconditioner
function z = eq_schur(params,x,dx)
    % Keep track of where the evaluation occurs
    persistent cache

    % Performance diagnostics
    global diagnostics

    % Here, we need to cache two elements due to the equality multiplier solve.
    % Basically, the equality multiplier solve during globalization requires a
    % solve at a new iterate.  If globalization accepts this point, we can
    % reuse this factorization.  However, if globalization rejects this point,
    % we want to use our old factorization.

    % Figure out if we match a cached element
    [cache iscached]=cache_search(cache,x);

    % If we don't have a match, cache a new factorization 
    if ~iscached 
        % Save the current location 
        cache{1}.x = x;

        % Exact Schur preconditioner
        if params.approx_schur==0
            % Factorize the total derivative of g'
            [q cache{1}.r] = qr(deriv(params,x)',0); 

        % Approximate Schur preconditioner
        else
            % Factorize the differential operator
            [cache{1}.l cache{1}.u cache{1}.p cache{1}.q cache{1}.r] = ...
                lu(op(params,x),'vector');
        end
        
        % Keep track that we did a new factorization 
        diagnostics.factorization_cached = diagnostics.factorization_cached+1;
    end
    
    % Solve the linear system 
    if params.approx_schur==0
        z = cache{1}.r\(cache{1}.r'\dx);
    else
        % Forward
        z=zeros(params.nx,1);
        z(cache{1}.q) = cache{1}.u\(cache{1}.l\(cache{1}.r(:,cache{1}.p)\dx));

        % Adjoint
        z = cache{1}.r(:,cache{1}.p)'\(cache{1}.l'\(cache{1}.u'\z(cache{1}.q)));
    end
end

% Adjust the forcing function with the pieces for the boundary conditions
function z = rhs(params,x)
    z = params.f - x(params.idx.k(1))*params.Ahat ...
        - x(params.idx.k(2))*params.Bhat;
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
    z = x(params.idx.k(1))*params.A+x(params.idx.k(2))*params.B;
end

% Grabs the derivative of the differential operator 
function z = op_p(which,params,x)
    if which==1
        z = params.A;
    else
        z = params.B;
    end
end

% Finds the total derivative of the equality constraints
function D = deriv(params,x)
    % Keep track of where the evaluation occurs
    persistent cache

    % Performance diagnostics
    global diagnostics
    
    % Figure out if we match a cached element
    [cache iscached]=cache_search(cache,x);

    % If we don't have a match, cache a new factorization 
    if ~iscached 
        % Save the current location 
        cache{1}.x = x;

        % Find the total derivative 
        cache{1}.D = [ ...
            op_p(1,params,x)*x(params.idx.u)-rhs_p(1,params,x) ...
            op_p(2,params,x)*x(params.idx.u)-rhs_p(2,params,x) ...
            op(params,x)];
        
        % Keep track that we cached a derivative 
        diagnostics.first_derivative_cached = ...
            diagnostics.first_derivative_cached+1;
    end

    % Return the derivative
    D = cache{1}.D;
end

% Finds the second total derivative adjoint of the equality constraints applied
% to the equality multiplier
function D2 = deriv2(params,x,dy)
    % Keep track of where the evaluation occurs
    persistent cache
    global diagnostics

    % Cache the total derivative when possible 
    if isempty(cache) || ~isequal(x,cache.x) || ~isequal(dy,cache.dy)
        % Save the current location 
        cache.x = x;
        cache.dy = dy;

        % Find the adjoint of the second derivative applied to the equality
        % multiplier
        cache.D2 = sparse(params.nx+2,params.nx+2);
        cache.D2(params.idx.k(1),params.idx.u) = dy'*op_p(1,params,x);
        cache.D2(params.idx.k(2),params.idx.u) = dy'*op_p(2,params,x);
        cache.D2(params.idx.u,params.idx.k(1)) = op_p(1,params,x)'*dy;
        cache.D2(params.idx.u,params.idx.k(2)) = op_p(2,params,x)'*dy;
        
        % Keep track that we cache a derivative 
        diagnostics.second_derivative_cached = ...
            diagnostics.second_derivative_cached+1;
    end

    % Return the derivative
    D2 = cache.D2;
end

% Prepares our cached element according to the following scheme
%
% 1.  Item not cached, copy first cached element to the second.  Return that no
%     cached item found.
%
% 2.  Item found in first cached element.  Return that cached item found.
% 
% 3.  Item found in second cached element.  Exchange first and second cached
%     elements.  Return that cached item found.
function [cache iscached] = cache_search(cache,x)
    % Determine what cached item matches x
    which = 0;
    if ~isempty(cache)
        for i=1:length(cache)
            if isequal(x,cache{i}.x)
                which = i;
                break;
            end
        end
    end

    % No items match
    if which==0
        iscached = 0; 
        if ~isempty(cache)
            cache{2} = cache{1};
        end

    % First item matches
    elseif which==1
        iscached = 1;

    % Second item matches
    elseif which==2
        iscached = 1;
        cache(2:-1:1)=cache;
    end
end

% Solves the discretized PDE without caching
function z = state_uncached(params,x,rhs)
    % Factorize the operator
    [l u p] = lu(op(params,x),'vector');
    
    % Solve the linear system 
    z = u\(l\rhs(p));
end
