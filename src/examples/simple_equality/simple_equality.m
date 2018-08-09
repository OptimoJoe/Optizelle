% Optimize a simple optimization problem with an optimal solution
% of (2-sqrt(2)/2,2-sqrt(2)/2).
function simple_equality(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('simple_equality <parameters>');
    end

    % Execute the optimization
    main(fname);
end

%---Objective0---
% Squares its input
function z = sq(x)
    z=x*x;
end

% Define a simple objective where
%
% f(x,y)=x^2+y^2
%
function self = MyObj()

    % Evaluation
    self.eval = @(x) sq(x(1))+sq(x(2));

    % Gradient
    self.grad = @(x) [ ...
        2.*x(1); ...
        2.*x(2)];

    % Hessian-vector product
    self.hessvec = @(x,dx) [ ...
        2.*dx(1); ...
        2.*dx(2)];
end
%---Objective1---

%---EqualityConstraint0---
% Define a simple equality constraint
%
% g(x,y)= [ (x-2)^2 + (y-2)^2 = 1 ]
%
function self = MyEq()

    % y=g(x)
    self.eval = @(x) [ ...
        sq(x(1)-2.)+sq(x(2)-2.)-1.];

    % y=g'(x)dx
    self.p = @(x,dx) [ ...
        2.*(x(1)-2.)*dx(1)+2.*(x(2)-2.)*dx(2)];

    % xhat=g'(x)*dy
    self.ps = @(x,dy) [ ...
        2.*(x(1)-2.)*dy(1); ...
        2.*(x(2)-2.)*dy(1)];

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) [ ...
        2.*dx(1)*dy(1); ...
        2.*dx(2)*dy(1) ];
end
%---EqualityConstraint1---

%---Preconditioner0---
% Define a Schur preconditioner for the equality constraints
function self = MyPrecon()
    self.eval=@(state,dy)dy(1)/sq(4.*(state.x(1)-2.)+4.*sq(state.x(2)-2.));
end
%---Preconditioner1---

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    %---State0---
    % Generate an initial guess
    x = [2.1;1.1];

    % Allocate memory for the equality multiplier
    y = [0.];

    % Create an optimization state
    state= Optizelle.EqualityConstrained.State.t(Optizelle.Rm,Optizelle.Rm,x,y);
    %---State1---

    %---Parameters0---
    % Read the parameters from file
    state = Optizelle.json.EqualityConstrained.read( ...
        Optizelle.Rm,Optizelle.Rm,fname,state);
    %---Parameters1---

    %---Functions0---
    % Create a bundle of functions
    fns=Optizelle.EqualityConstrained.Functions.t;
    fns.f=MyObj();
    fns.g=MyEq();
    fns.PSchur_left=MyPrecon();
    %---Functions1---

    %---Solver0---
    % Solve the optimization problem
    state = Optizelle.EqualityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);
    %---Solver1---

    %---Extract0---
    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));
    %---Extract1---

    % Write out the final answer to file
    Optizelle.json.EqualityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,'solution.json',state);
end
