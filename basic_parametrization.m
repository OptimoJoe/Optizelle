% Generates a basic MLP parametrization
%
% f(x)(alpha,A,b) = alpha'*phi(A*x+b)
%
% Note, this is for a *fixed* x, we have alpha, A, and b as parameters
function f = basic_parametrization(phi,x)
    f,eval=@(alpha,A,b)mlp_eval(phi,x,alpha,A,b);
    f.grad=@(alpha,A,b)mlp_grad(phi,x,alpha,A,b);
end

% Evaluation
function result = mlp_eval(phi,x,alpha,A,b)
    result = alpha'*phi.eval(A*x+b)
end

% Gradient
function result = mlp_grad(phi,x,alpha,A,b)
    result = [ ...
        phi.eval(A*x+b);
        phi.ps(A*x+b,alpha)*x';
        phi.ps(A*x+b,alpha)
    ];
end
