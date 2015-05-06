% Generates a the objective function 
%
% J(alpha,A,b) = .5 * norm(alpha'*phi(A*x+b) - y)
%
% where the sizes of each variable and parameter are 
%
% alpha : nhidden
% A : nhidden x ninput
% b : nhidden
% x : ninput x nsamples
% y : nsamples
function f = generate_objective(phi,lens,x,y)
    f.eval=@(xx)mlp_eval(phi,lens,x,y,xx);
    f.grad=@(xx)mlp_grad(phi,lens,x,y,xx);
    f.hessvec=@(xx,dxx)mlp_hessvec(phi,lens,x,y,xx,dxx);
end

% Evaluation
function result = mlp_eval(phi,lens,x,y,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Duplicate by rows
    bb = @(x)repmat(b,1,size(x,2));

    % Evaluate the objective
    result = 0.5 * norm(phi.eval(A*x+bb(x))'*alpha-y,2)^2;
end

% Gradient
function result = mlp_grad(phi,lens,x,y,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Duplicate by rows
    bb = @(x)repmat(b,1,size(x,2));

    % Vector of ones
    e = ones(size(x,2),1);

    % Create an empty gradient
    result = zeros(size(xx));

    % Find the gradient in the different directions
    result_alpha = ...
        phi.eval(A*x+bb(x))*(phi.eval(A*x+bb(x))'*alpha-y);

    result_A = ...
        phi.ps(A*x+bb(x),alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*x';

    result_b = ...
        phi.ps(A*x+bb(x),alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*e;

    % Find the gradient 
    result = ...
        lens.alpha.set(result_alpha, ...
        lens.A.set(result_A, ...
        lens.b.set(result_b, ...
        result)));
        
end

% Hessian-vector product 
function result = mlp_hessvec(phi,lens,x,y,xx,dxx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);
    
    dalpha = lens.alpha.get(dxx);
    dA = lens.A.get(dxx);
    db = lens.b.get(dxx);

    % Duplicate by rows
    bb = @(x)repmat(b,1,size(x,2));
    dbb = @(x)repmat(db,1,size(x,2));

    % Vector of ones
    e = ones(size(x,2),1);

    % Create an empty Hessian-vector product 
    result = zeros(size(xx));

    % Find the curvature in the different directions 
    result_alpha_alpha = phi.eval(A*x+bb(x))*(phi.eval(A*x+bb(x))'*dalpha);

    result_alpha_A = ... 
            phi.p(A*x+bb(x),dA*x)*(phi.eval(A*x+bb(x))'*alpha-y) +...
            phi.eval(A*x+bb(x))*(phi.p(A*x+bb(x),dA*x)'*alpha);

    result_alpha_b = ... 
            phi.p(A*x+bb(x),dbb(x))*(phi.eval(A*x+bb(x))'*alpha-y) +...
            phi.eval(A*x+bb(x))*(phi.p(A*x+bb(x),dbb(x))'*alpha);

    result_A_alpha = ...
            phi.ps(A*x+bb(x),dalpha*(phi.eval(A*x+bb(x))'*alpha-y)' + ...
                             alpha*(phi.eval(A*x+bb(x))'*dalpha)')*x';
    
    result_A_A = ...
        phi.pps(A*x+bb(x),dA*x,alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*x' + ...
        phi.ps(A*x+bb(x),alpha*(phi.p(A*x+bb(x),dA*x)'*alpha)')*x';
    
    result_A_b = ...
        phi.pps(A*x+bb(x),dbb(x),alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*x' + ...
        phi.ps(A*x+bb(x),alpha*(phi.p(A*x+bb(x),dbb(x))'*alpha)')*x';
    
    result_b_alpha = ...
        phi.ps(A*x+bb(x),dalpha*(phi.eval(A*x+bb(x))'*alpha-y)' + ...
                         alpha*(phi.eval(A*x+bb(x))'*dalpha)')*e;
    
    result_b_A = ...
            phi.pps(A*x+bb(x),dA*x,alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*e +...
            phi.ps(A*x+bb(x),alpha*(phi.p(A*x+bb(x),dA*x)'*alpha)')*e;
    
    result_b_b = ...
            phi.pps(A*x+bb(x),dbb(x),alpha*(phi.eval(A*x+bb(x))'*alpha-y)')*e...
            + ...
            phi.ps(A*x+bb(x),alpha*(phi.p(A*x+bb(x),dbb(x))'*alpha)')*e;

    % Find the Hessian-vector product
    result = ...
        lens.alpha.set(result_alpha_alpha + result_alpha_A + result_alpha_b, ...
        lens.A.set(result_A_alpha + result_A_A + result_A_b, ...
        lens.b.set(result_b_alpha + result_b_A + result_b_b, ...
        result)));
end
