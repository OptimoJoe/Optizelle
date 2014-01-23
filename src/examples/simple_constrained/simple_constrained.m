% Optimize a simple optimization problem with an optimal solution
% of (1/3,1/3)
function simple_constrained(fname)
    % Read in the name for the input file
    if nargin ~=2
        error('simple_constrained <parameters>');
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
function self = MyObj(self)

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

% Define a simple equality
%
% g(x,y)= [ x + 2y = 1 ] 
%
function self = MyEq(self)

    % y=g(x) 
    self.eval = @(x) [x(1)+2.*x(2)-1.];

    % y=g'(x)dx
    self.p = @(x,dx) [dx(1)+2.*dx(2)];

    % z=g'(x)*dy
    self.ps = @(x,dy) [
        dy(1);
        2.*dy(1)];

    % z=(g''(x)dx)*dy
    self.pps = @(x,dy,dy) zeros(2,1); 
end

% Define simple inequalities 
%
% h(x,y)= [ 2x + y >= 1 ] 
%
function self = MyIneq(self)

    % y=h(x) 
    self.eval = @(x) [
        2.*x(1)+x(2)-1];

    % y=h'(x)dx
    self.p = @(x,dx) [
        2.*dx(1)+dx(2)];

    % z=h'(x)*dy
    self.ps = @(x,dy) [
        2.*dy(1)
        dy(1)];

    % z=(h''(x)dx)*dy
    self.pps = @(x,dy,dy) [ 0. ]; 
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    Optizelle = setupOptizelle ();

    % Generate an initial guess 
    x = [2.1;1.1];

    % Allocate memory for the equality multiplier 
    y = [0.];

    % Allocate memory for the inequality multiplier 
    z = [0.];

    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,x,y,z);

    % Read the parameters from file
    Optizelle.json.Constrained.read( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fname,state);

    % Create a bundle of functions
    fns=Optizelle.Constrained.Functions.t;
    fns.f=MyObj(Optizelle.ScalarValuedFunction);
    fns.g=MyEq(Optizelle.ScalarValuedFunction);
    fns.h=MyIneq(Optizelle.ScalarValuedFunction);

    % Solve the optimization problem
    Optizelle.Constrained.Algorithms.getMin(
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));

    % Write out the final answer to file
    Optizelle.json.Constrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm, ...
        Optizelle.Messaging(),'solution.json',state);
end
