% Generate an interpolatory function
%
% f(x) = alpha'*phi(A*x)
%
function f=generate_interpolant(phi,alpha,A,b)
    f.eval = @(x)alpha'*phi(A*x+b);
end
