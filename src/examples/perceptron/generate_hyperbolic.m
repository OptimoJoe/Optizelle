% Creates a vectorized version of the hyperbolic tangent function along
% with its derivatives.
function f = generate_hyperbolic()
    f.eval = @(x)tanh_eval(x);
    f.p = @(x,dx)tanh_p(x).*dx;
    f.ps = @(x,dy)tanh_p(x).*dy;
    f.pps = @(x,dx,dy)tanh_p2(x).*dx.*dy;
end

% Evaluation
function result = tanh_eval(x)
    result = tanh(x);
end

% Derivative 
function result = tanh_p(x) 
    result = sech(x).^2;
end

% Second derivative
function result = tanh_p2(x) 
    result = -2.0*sech(x).^2.*tanh(x);
end
