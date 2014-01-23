% Optimize a simple optimization problem with an optimal solution 
% of (2-sqrt(2)/2,2-sqrt(2)/2).
function simple_equality(fname)
    % Read in the name for the input file
    if nargin ~=2
        error('simple_equality <parameters>');
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
% f(x,y)=x^2+y^2
%
function self = MyObj(self)

    % Evaluation 
    self.eval = @(x) sq(x(1))+sq(x(2));

    % Gradient
    self.grad = @(x) [
        2.*x(1);
        2.*x(2)];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
        2.*dx(1);
        2.*dx(2)];

% Define a simple equality constraint
%
% g(x,y)= [ (x-2)^2 + (y-2)^2 = 1 ] 
%
function self = MyEq(self)

    % y=g(x) 
    self.eval = @(x) [
        sq(x(1)-2.)+sq(x(2)-2.)-1.];

    % y=g'(x)dx
    self.p = @(x,dx) [
        y(1) = 2.*(x(1)-2.)*dx(1)+2.*(x(2)-2.)*dx(2)];

    % z=g'(x)*dy
    self.ps = @(x,dy) [
        z(1) = 2.*(x(1)-2.)*dy(1)
        z(2) = 2.*(x(2)-2.)*dy(1)];

    % z=(g''(x)dx)*dy
    self.pps = @(x,dy,dy) [
        2.*dx(1)*dy(1);
        2.*dx(2)*dy(1) ];

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    Optizelle = setupOptizelle ();

    % Generate an initial guess 
    x = [2.1;1.1];

    % Allocate memory for the equality multiplier 
    y = [0.];

    % Create an optimization state
    state=Optizelle.EqualityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,x,y);

    % Read the parameters from file
    Optizelle.json.EqualityConstrained.read(
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fname,state);

    % Create a bundle of functions
    fns=Optizelle.EqualityConstrained.Functions.t;
    fns.f=MyObj(Optizelle.ScalarValuedFunction);
    fns.g=MyEq(Optizelle.VectorValuedFunction);

    % Solve the optimization problem
    Optizelle.EqualityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));

    % Write out the final answer to file
    Optizelle.json.EqualityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,'solution.json',state);
end
