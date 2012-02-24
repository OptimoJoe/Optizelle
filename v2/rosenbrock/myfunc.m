% Rosenbrock function
function z=myfunc(x)
    z=(1-x(1))^2+100*(x(2)-x(1)^2)^2;
