% Hessian vector product for the Rosenbrock function
function [z H]=myhess(x,eta)
    H=[1200*x(1)^2-400*x(2)+2 -400*x(1);-400*x(1) 200];
    z=H*eta;
