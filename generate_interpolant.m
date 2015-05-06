% Generate an interpolatory function
%
% f(x) = alpha'*phi(A*x+b)
%
% This works with vectorized inputs where size(x)=[ninput,nsamples].  In this
% case,
%
% f(x) : 1 x nsamples
% grad f(x) : ninput x nsamples
% hess f(x) dx : ninput x nsamples
%
function f=generate_interpolant(phi,alpha,A,b)
    % Figure out the size of the problem
    [nhidden ninput] = size(A);

    % Duplicate b by rows.  This is necessary for the vectorization.
    bb = @(x)repmat(b,1,size(x,2));
    aalpha = @(x)repmat(alpha,1,size(x,2));

    % Create the interpolant
    f.eval = @(x)alpha'*phi.eval(A*x+bb(x));
    f.grad = @(x)A'*phi.ps(A*x+bb(x),aalpha(x));
    f.hessvec = @(x,dx)A'*phi.pps(A*x+bb(x),A*dx,aalpha(x));
end
