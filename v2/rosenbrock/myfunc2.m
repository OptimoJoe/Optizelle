% Rosenbrock function
function z=myfunc(x,y)
    z=(1-x).^2+100*(y-x.^2).^2;
