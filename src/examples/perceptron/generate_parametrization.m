% Generates a parametrization
%
% f(x)(alpha,A,b) = alpha'*phi(A*x+b)+beta
%
% This works with vectorized inputs where size(x)=[ninput,nsamples].  In this
% case,
%
% f(alpha,A,b) : 1 x nsamples
% grad f(alpha,A,b) : (nhidden + 1 + nhidden*ninput + nhidden) x nsamples
% hess f(alpha,A,b) dx : (nhidden + 1 + nhidden*ninput + nhidden) x nsamples
%
function [f scaling] = generate_parametrization(phi,lens,x,scaling)
    % Generate the scaling operators
    scaling_x = generate_scaling_operator( ...
        scaling.x.from.min,scaling.x.from.max, ...
        scaling.x.to.min,scaling.x.to.max);

    scaling_y_inv = generate_scaling_operator( ...
        scaling.y.to.min,scaling.y.to.max, ...
        scaling.y.from.min,scaling.y.from.max);

    % Scale the inputs
    x = scaling_x.eval(x);

    f.eval=@(xx)scaling_y_inv.eval(mlp_eval(phi,lens,x,xx));
    f.grad=@(xx)scaling_y_inv.ps(1,1)*mlp_grad(phi,lens,x,xx);
    f.hessvec=@(xx,dxx)scaling_y_inv.ps(1,1)*mlp_hessvec(phi,lens,x,xx,dxx);
end

% Evaluation
function result = mlp_eval(phi,lens,x,xx)
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
    bbeta = repmat(beta,1,nsamples);

    % Evaluate the MLP
    result = alpha'*phi.eval(A*x+bb)+bbeta;
end

% Gradient
function result = mlp_grad(phi,lens,x,xx)
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
    aalpha = repmat(alpha,1,nsamples);

    % Calculate the gradient
    result = [ ...
        phi.eval(A*x+bb);
        ones(1,nsamples);
        myouter(phi.ps(A*x+bb,aalpha),x);
        phi.ps(A*x+bb,aalpha)
    ];
end

% Hessian-vector product
function result = mlp_hessvec(phi,lens,x,xx,dxx)
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
    aalpha = repmat(alpha,1,nsamples);

    dbb = repmat(db,1,nsamples);
    daalpha = repmat(dalpha,1,nsamples);

    % Calculate the Hessian-vector product
    hv_alpha_alpha = zeros(nhidden,nsamples);
    hv_alpha_beta = zeros(nhidden,nsamples);
    hv_alpha_A = phi.p(A*x+bb,dA*x);
    hv_alpha_b = phi.p(A*x+bb,dbb);

    hv_A_alpha = myouter(phi.ps(A*x+bb,daalpha),x);
    hv_A_beta = zeros(ninput*nhidden,nsamples);
    hv_A_A = myouter(phi.pps(A*x+bb,dA*x,aalpha),x);
    hv_A_b = myouter(phi.pps(A*x+bb,dbb,aalpha),x);

    hv_b_alpha = phi.ps(A*x+bb,daalpha);
    hv_b_beta = zeros(nhidden,nsamples);
    hv_b_A = phi.pps(A*x+bb,dA*x,aalpha);
    hv_b_b = phi.pps(A*x+bb,dbb,aalpha);

    result = [ ...
        hv_alpha_alpha + hv_alpha_beta + hv_alpha_A + hv_alpha_b;
        zeros(1,nsamples);
        hv_A_alpha + hv_A_beta + hv_A_A + hv_A_b;
        hv_b_alpha + hv_b_beta + hv_b_A + hv_b_b
    ];
end

% Takes the outer product between the corresponding columns and stacks the
% result
function result = myouter(A,B)
    % Figure out the sizes of our matrices.  We assume they have the same
    % number of columns.
    [m n] = size(A);
    p = size(B,1);

    % Allocate memory for the result
    result=zeros(m*p,n);

    % Loop over the columns of A
    for j=1:n
        outer = A(:,j)*B(:,j)';
        result(:,j) = outer(:);
    end
end
