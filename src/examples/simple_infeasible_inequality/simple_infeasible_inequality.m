% Optimize a simple optimization problem with an optimal solution of (1/3,1/3),
%
% min x + y
% st  x + 2y >= 1
%     2x + y >= 1
%
% Now, in the case we don't have a starting feasible solution, we can play
% a reformulation trick that adds two scalar variables and allows us to find
% a strictly feasible solution.  Namely,
%
% Note, most of the time, we're much better off just adding slack variables.
% Basically, this trick is only worthwhile when we don't have a linear system
% solver for the equality constraints added from the slacks since this method
% only adds a single equality constraint.

% min x + y
% st  x + 2y >= 1 - z
%     2x + y >= 1 - z
%     epsilon >= w
%     z = w
function simple_infeasible_inequality(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('simple_infeasible_inequality <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Define an objective where
%
% f(x,y,z,w)=x+y
%
function self = MyObj()

    % Evaluation
    self.eval = @(x) x(1)+x(2);

    % Gradient
    self.grad = @(x) [ ...
        1.; ...
        1.; ...
        0.; ...
        0.];

    % Hessian-vector product
    self.hessvec = @(x,dx) zeros(4,1);
end

% Define a single equality where
%
% g(x,y,z,w) = z - w = 0
%
function self = MyEq()

    % y=g(x)
    self.eval = @(x) [x(3)-x(4)];

    % y=g'(x)dx
    self.p = @(x,dx) [dx(3)-dx(4)];

    % xhat=g'(x)*dy
    self.ps = @(x,dy) [ ...
        0.; ...
        0.; ...
        dy(1); ...
        -dy(1)];

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) zeros(4,1);
end

% Define some inequalities where
%
% h(x,y,z,w) = [ x + 2y >= 1 - z  ]
%              [ 2x + y >= 1 - z  ]
%              [ epsilon >= w     ]
%
function self = MyIneq(epsilon)

    % z=h(x)
    self.eval = @(x) [ ...
        x(1)+2.*x(2)+x(3)-1.; ...
        2.*x(1)+x(2)+x(3)-1.; ...
        epsilon-x(4)];

    % z=h'(x)dx
    self.p = @(x,dx) [ ...
        dx(1)+2.*dx(2)+dx(3); ...
        2.*dx(1)+dx(2)+dx(3); ...
        -dx(4)];

    % xhat=h'(x)*dz
    self.ps = @(x,dz) [ ...
        dz(1)+2.*dz(2); ...
        2.*dz(1)+dz(2); ...
        dz(1)+dz(2); ...
        -dz(3)];

    % xhat=(h''(x)dx)*dz
    self.pps = @(x,dx,dy) zeros(4,1);
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Set the amount of infeasibility that we want to allow
    epsilon = 1e-8;

    % Generate an initial guess for the primal
    x = [0.;0.;5.;-5.];

    % Generate a vector for the equality multiplier
    y = [0.];

    % Generate a vector for the inequality multiplier
    z = [0.;0.;0.];

    % Create an optimization state
    state=Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z);

    % Read the parameters from file
    state=Optizelle.json.Constrained.read( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,fname,state);

    % Create a bundle of functions
    fns = Optizelle.Constrained.Functions.t;
    fns.f = MyObj();
    fns.g = MyEq();
    fns.h = MyIneq(epsilon);

    % Solve the optimization problem
    state = Optizelle.Constrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
        fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));

    % Write out the final answer to file
    Optizelle.json.Constrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,'solution.json',state);
end
