% Generates a the objective function 
%
% J(alpha,beta,A,b) = .5 * norm(phi(A*x+b*e')' * alpha + beta*e - y')
%
% where the sizes of each variable and parameter are 
%
% alpha : nhidden
% beta : 1
% A : nhidden x ninput
% b : nhidden
% x : ninput x nsamples
% y : 1 x nsamples
function [f scaling] = generate_objective(phi,lens,x,y,scaling)
    % Generate the scaling operators 
    scaling_x = generate_scaling_operator( ...
        scaling.x.from.min,scaling.x.from.max, ...
        scaling.x.to.min,scaling.x.to.max);

    scaling_y = generate_scaling_operator( ...
        scaling.y.from.min,scaling.y.from.max, ...
        scaling.y.to.min,scaling.y.to.max);

    % Scale the inputs and the targets 
    x = scaling_x.eval(x);
    y = scaling_y.eval(y);

    % Find the objective
    f.eval=@(xx)mlp_eval(phi,lens,x,y,xx);
    f.grad=@(xx)mlp_grad(phi,lens,x,y,xx);
    f.hessvec=@(xx,dxx)mlp_hessvec(phi,lens,x,y,xx,dxx);
end

% Evaluation
function result = mlp_eval(phi,lens,x,y,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    beta = lens.beta.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Grab sizes
    [ninput nsamples]=size(x);
    nhidden=length(b);

    % Expand terms 
    bb = repmat(b,1,nsamples);
    bbeta = repmat(beta,nsamples,1);

    % Evaluate the objective
    result = 0.5 * norm(phi.eval(A*x+bb)'*alpha+bbeta-y',2)^2;
end

% Gradient
function result = mlp_grad(phi,lens,x,y,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    beta = lens.beta.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Grab sizes
    [ninput nsamples]=size(x);
    nhidden=length(b);

    % Expand terms 
    bb = repmat(b,1,nsamples);
    bbeta = repmat(beta,nsamples,1);

    % Vector of ones
    e = ones(nsamples,1);

    % Create an empty gradient
    result = zeros(size(xx));

    % Find the gradient in the different directions
    result_alpha = ...
        phi.eval(A*x+bb)*(phi.eval(A*x+bb)'*alpha+bbeta-y');
    
    result_beta = ...
        sum(phi.eval(A*x+bb)'*alpha+bbeta-y');
    
    result_A = ...
        phi.ps(A*x+bb,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')*x';

    result_b = ...
        phi.ps(A*x+bb,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')*e;

    % Find the gradient 
    result = ...
        lens.alpha.set(result_alpha, ...
        lens.beta.set(result_beta, ...
        lens.A.set(result_A, ...
        lens.b.set(result_b, ...
        result))));
end

% Hessian-vector product 
function result = mlp_hessvec(phi,lens,x,y,xx,dxx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    beta = lens.beta.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);
    
    dalpha = lens.alpha.get(dxx);
    dbeta = lens.beta.get(dxx);
    dA = lens.A.get(dxx);
    db = lens.b.get(dxx);

    % Grab sizes
    [ninput nsamples]=size(x);
    nhidden=length(b);

    % Expand terms 
    bb = repmat(b,1,nsamples);
    bbeta = repmat(beta,nsamples,1);
    
    dbb = repmat(db,1,nsamples);
    dbbeta = repmat(dbeta,nsamples,1);

    % Vector of ones
    e = ones(nsamples,1);

    % Create an empty Hessian-vector product 
    result = zeros(size(xx));

    % Find the curvature in the different directions 

    % alpha
    result_alpha_alpha = phi.eval(A*x+bb)*(phi.eval(A*x+bb)'*dalpha);
    
    result_alpha_beta = phi.eval(A*x+bb)*(e*dbeta);

    result_alpha_A = ... 
            phi.p(A*x+bb,dA*x)*(phi.eval(A*x+bb)'*alpha+bbeta-y') +...
            phi.eval(A*x+bb)*(phi.p(A*x+bb,dA*x)'*alpha);

    result_alpha_b = ... 
            phi.p(A*x+bb,dbb)*(phi.eval(A*x+bb)'*alpha+bbeta-y') +...
            phi.eval(A*x+bb)*(phi.p(A*x+bb,dbb)'*alpha);
    
    result_alpha = result_alpha_alpha + result_alpha_beta + result_alpha_A ...
        + result_alpha_b;
   
    % beta
    result_beta_alpha = ...
        sum(phi.eval(A*x+bb)'*dalpha);
    
    result_beta_beta = sum(dbbeta); 
    
    result_beta_A = sum(phi.p(A*x+bb,dA*x)'*alpha);
    
    result_beta_b = sum(phi.p(A*x+bb,dbb)'*alpha);
    
    result_beta = result_beta_alpha + result_beta_beta + result_beta_A ...
        + result_beta_b;

    % A
    result_A_alpha = ...
        phi.ps(A*x+bb,dalpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')' + ...
                         alpha*(phi.eval(A*x+bb)'*dalpha)')*x';

    result_A_beta = ...
        phi.ps(A*x+bb,alpha*dbbeta')*x';
    
    result_A_A = ...
        phi.pps(A*x+bb,dA*x,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')*x' + ...
        phi.ps(A*x+bb,alpha*(phi.p(A*x+bb,dA*x)'*alpha)')*x';
    
    result_A_b = ...
        phi.pps(A*x+bb,dbb,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')*x' +...
        phi.ps(A*x+bb,alpha*(phi.p(A*x+bb,dbb)'*alpha)')*x';

    result_A = result_A_alpha + result_A_beta + result_A_A + result_A_b;
    
    % b
    result_b_alpha = ...
        phi.ps(A*x+bb,dalpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')' + ...
                         alpha*(phi.eval(A*x+bb)'*dalpha)')*e;
    
    result_b_beta = ...
        phi.ps(A*x+bb,alpha*dbbeta')*e;
    
    result_b_A = ...
            phi.pps(A*x+bb,dA*x,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')*e+...
            phi.ps(A*x+bb,alpha*(phi.p(A*x+bb,dA*x)'*alpha)')*e;
    
    result_b_b = ...
            phi.pps(A*x+bb,dbb,alpha*(phi.eval(A*x+bb)'*alpha+bbeta-y')')...
                *e+ ...
            phi.ps(A*x+bb,alpha*(phi.p(A*x+bb,dbb)'*alpha)')*e;

    result_b = result_b_alpha + result_b_beta + result_b_A + result_b_b;

    % Find the Hessian-vector product
    result = ...
        lens.alpha.set(result_alpha, ...
        lens.beta.set(result_beta, ...
        lens.A.set(result_A, ...
        lens.b.set(result_b, ...
        result))));
end
