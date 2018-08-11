% Creates a vectorized version of the logistic function
%
% f = L / (1 + exp(-k(x-x0)))
%
% along with its derivatives
function f = generate_logistic(L,k,x0)
    f.eval = @(x)logistic_eval(L,k,x0,x);
    f.p = @(x,dx)logistic_p(L,k,x0,x).*dx;
    f.ps = @(x,dy)logistic_p(L,k,x0,x).*dy;
    f.pps = @(x,dx,dy)logistic_p2(L,k,x0,x).*dx.*dy;
end

% Evaluation
function result = logistic_eval(L,k,x0,x)
    result = L./(1+exp(-k*(x-x0)));
end

% Derivative
function result = logistic_p(L,k,x0,x)
    result = L*k*exp(-k*(x-x0))./(exp(-k*(x-x0))+1).^2;
end

% Second derivative
function result = logistic_p2(L,k,x0,x)
    result= -L*k^2*(exp(k*x)-exp(k*x0)).*exp(k*(x+x0))./(exp(k*x)+exp(k*x0)).^3;
end
