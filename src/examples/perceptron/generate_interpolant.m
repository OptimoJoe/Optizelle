% Generate an interpolatory function
%
% f(alpha,A,b)(x) = alpha'*phi(A*x+b)+beta
%
% This works with vectorized inputs where size(x)=[ninput,nsamples].  In this
% case,
%
% f(x) : 1 x nsamples
% grad f(x) : ninput x nsamples
% hess f(x) dx : ninput x nsamples
%
function f=generate_interpolant(phi,lens,xx,scaling)
    % Generate the scaling operators
    scaling_x = generate_scaling_operator( ...
        scaling.x.from.min,scaling.x.from.max, ...
        scaling.x.to.min,scaling.x.to.max);

    scaling_y_inv = generate_scaling_operator( ...
        scaling.y.to.min,scaling.y.to.max, ...
        scaling.y.from.min,scaling.y.from.max);

    % Grab our parameters
    alpha = lens.alpha.get(xx);
    beta = lens.beta.get(xx);
    A = lens.A.get(xx);
    b = lens.b.get(xx);

    % Expand terms
    bb = @(x)repmat(b,1,size(x,2));
    aalpha = @(x)repmat(alpha,1,size(x,2));
    bbeta = @(x)repmat(beta,1,size(x,2));

    % Create the interpolant
    f.eval = @(x)scaling_y_inv.eval( ...
        alpha'*phi.eval(A*scaling_x.eval(x)+bb(x))+bbeta(x));
    f.grad = @(x)scaling_y_inv.ps(1,1) * ...
        scaling_x.ps(x,A'*phi.ps(A*scaling_x.eval(x)+bb(x),aalpha(x)));
    f.hessvec = @(x,dx)scaling_y_inv.ps(1,1) * scaling_x.ps( ...
        x,A'*phi.pps(A*scaling_x.eval(x)+bb(x),A*scaling_x.p(x,dx),aalpha(x)));
end
