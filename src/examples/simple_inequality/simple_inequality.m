% Optimize a simple optimization problem with an optimal solution
% of (1/3,1/3)
function simple_inequality(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('simple_inequality <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Define a simple objective where 
% 
% f(x,y)=(x+1)^2+(y+1)^2
%
function self = MyObj()

    % Evaluation 
    self.eval = @(x) sq(x(1)+1.)+sq(x(2)+1.);

    % Gradient
    self.grad = @(x) [
        2.*x(1)+2.;
        2.*x(2)+2.];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
        2.*dx(1);
        2.*dx(2)];
end

% Define simple inequalities 
%
% h(x,y)= [ x + 2y >= 1 ] 
%         [ 2x + y >= 1 ] 
%
function self = MyIneq()

    % z=h(x) 
    self.eval = @(x) [
        x(1)+2.*x(2)-1. ;
        2.*x(1)+x(2)-1. ];

    % z=h'(x)dx
    self.p = @(x,dx) [
        dx(1)+2.*dx(2) ;
        2.*dx(1)+dx(2) ];

    % xhat=h'(x)*dz
    self.ps = @(x,dz) [
        dz(1)+2.*dz(2) ;
        2.*dz(1)+dz(2) ];

    % xhat=(h''(x)dx)*dz
    self.pps = @(x,dx,dz) zeros(2,1); 
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess 
    x = [2.1;1.1];

    % Allocate memory for the inequality multiplier 
    z = [0.;0.];

    % Create an optimization state
    state=Optizelle.InequalityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,x,z);

    % Read the parameters from file
    state=Optizelle.json.InequalityConstrained.read( ...
        Optizelle.Rm,Optizelle.Rm,fname,state);

    % Create a bundle of functions
    fns=Optizelle.InequalityConstrained.Functions.t;
    fns.f=MyObj();
    fns.h=MyIneq();

    % Solve the optimization problem
    state=Optizelle.InequalityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));
    
    % Write out the final answer to file
    Optizelle.json.InequalityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,'solution.json',state);
end
