% Generates a parametrization 
%
% f(x)(alpha,A,b) = alpha'*phi(A*x+b)
%
% This works with vectorized inputs where size(x)=[ninput,nsamples].  In this
% case,
%
% f(alpha,A,b) : 1 x nsamples
% grad f(alpha,A,b) : (nhidden + nhidden*ninput + nhidden) x nsamples
% hess f(alpha,A,b) dx : (nhidden + nhidden*ninput + nhidden) x nsamples
%
function f = generate_parametrization(phi,lens,x)
    f.eval=@(xx)mlp_eval(phi,lens,x,xx);
    f.grad=@(xx)mlp_grad(phi,lens,x,xx);
    f.hessvec=@(xx,dxx)mlp_hessvec(phi,lens,x,xx,dxx);
end

% Evaluation
function result = mlp_eval(phi,lens,x,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);
    
    % Duplicate by rows
    bb = @(x)repmat(b,1,size(x,2));

    % Evaluate the MLP
    result = alpha'*phi.eval(A*x+bb(x));
end

% Gradient
function result = mlp_grad(phi,lens,x,xx)
    % Grab our parameters 
    alpha = lens.alpha.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Duplicate by rows
    bb = @(x)repmat(b,1,size(x,2));
    aalpha = @(x)repmat(alpha,1,size(x,2));

    % Calculate the gradient
    result = [ ...
        phi.eval(A*x+bb(x));
        myouter(phi.ps(A*x+bb(x),aalpha(x)),x);
        phi.ps(A*x+bb(x),aalpha(x))
    ];
end

% Hessian-vector product 
function result = mlp_hessvec(phi,lens,x,xx,dxx)
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
    aalpha = @(x)repmat(alpha,1,size(x,2));
    daalpha = @(x)repmat(dalpha,1,size(x,2));

    % Calculate the Hessian-vector product
    hv_alpha_alpha = zeros(size(alpha,1),size(x,2));
    hv_alpha_A = phi.p(A*x+bb(x),dA*x);
    hv_alpha_b = phi.p(A*x+bb(x),dbb(x));

    hv_A_alpha = myouter(phi.ps(A*x+bb(x),daalpha(x)),x);
    hv_A_A = myouter(phi.pps(A*x+bb(x),dA*x,aalpha(x)),x);
    hv_A_b = myouter(phi.pps(A*x+bb(x),dbb(x),aalpha(x)),x);
    
    hv_b_alpha = phi.ps(A*x+bb(x),daalpha(x));
    hv_b_A = phi.pps(A*x+bb(x),dA*x,aalpha(x));
    hv_b_b = phi.pps(A*x+bb(x),dbb(x),aalpha(x));

    result = [ ...
        hv_alpha_alpha + hv_alpha_A + hv_alpha_b;
        hv_A_alpha + hv_A_A + hv_A_b;
        hv_b_alpha + hv_b_A + hv_b_b
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
