% A problem that helps us determine how to scale inequality constrained
% optimization problems
function inequality_scaling(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('inequality_scaling <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Define a simple objective where
%
% f(x) = 0.5 || x - c ||^2
%
function self = MyObj(c)

    % Evaluation
    self.eval = @(x) 0.5*norm(x-c,2)^2;

    % Gradient
    self.grad = @(x) x-c;

    % Hessian-vector product
    self.hessvec = @(x,dx) dx;
end

% Define simple inequalities
%
% h(x) = x - lb
%
function self = MyIneq(lb)

    % z=h(x)
    self.eval = @(x) x-lb;

    % z=h'(x)dx
    self.p = @(x,dx) dx;

    % xhat=h'(x)*dz
    self.ps = @(x,dz) dz;

    % xhat=(h''(x)dx)*dz
    self.pps = @(x,dx,dz) zeros(size(x),1);
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Set the size
    m = 10;

    % Generate an initial guess
    x = ones(m,1)+10.^(-(1:m))';

    % Allocate memory for the inequality multiplier
    z = zeros(m,1);

    % Create the center of the objective function
    c = -ones(m,1);

    % Create the lower bound for the problem
    lb = ones(m,1);

    % Create an optimization state
    state=Optizelle.InequalityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,x,z);

    % Read the parameters from file
    state=Optizelle.json.InequalityConstrained.read( ...
        Optizelle.Rm,Optizelle.Rm,fname,state);

    % Create a bundle of functions
    fns=Optizelle.InequalityConstrained.Functions.t;
    fns.f=MyObj(c);
    fns.h=MyIneq(lb);

    % Solve the optimization problem
    state=Optizelle.InequalityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: [\n');
    for i=1:m
        fprintf('%1.16e\n',state.x(i));
    end
    fprintf(']\n');

    % Write out the final answer to file
    Optizelle.json.InequalityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,'solution.json',state);
end
