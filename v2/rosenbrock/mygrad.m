% Gradient for the Rosenbrock function
function z=mygrad(x)
    z=[-400*x(1)*(x(2)-x(1)^2) - 2 * (1 - x(1));200*(x(2)-x(1)^2)];
