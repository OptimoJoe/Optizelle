% Solve the problem %
% min x + y + z
% st  [ x y; y z] >= 0 <==> xz >= y^2
%     y >=Q z <==> y >= |z|
%     z >= 1
%
% which should have an optimal solution of (1,1,1).  Basically, we have
% a problem with a semidefinite, quadratic, and linear cone.  Then, we're
% going to run a few iterations, stop the optimization, write the restart,
% load the restart, and then continue with the optimization.  This is basically
% a way to both show a bunch of features and add another unit test.

function sql_restart(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('sql_restart <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Define a simple objective where 
% 
% f(x,y,z)=x+y+z
%
function self = MyObj()

    % Evaluation 
    self.eval = @(x) sum(x);

    % Gradient
    self.grad = @(x) ones(size(x)); 

    % Hessian-vector product
    self.hessvec = @(x,dx) zeros(size(x)); 
end

% Define a simple SQL inequality 
%
% h(x,y,z) = [ x y ] >=S 0
%            [ y z ]
%            y >=Q z
%            z >=L 1
%
function self = MyIneq()

    % z=h(x) 
    self.eval = @(x)MyIneq_eval(x);

    % z=h'(x)dx
    self.p = @(x,dx)MyIneq_p(x,dx);

    % xhat=h'(x)*dz
    self.ps = @(x,dz)MyIneq_ps(x,dz);

    % xhat=(h''(x)dx)*dz
    self.pps = @(x,dx,dz) zeros(size(x)); 
end

% z=h(x) 
function z=MyIneq_eval(x)
    global Optizelle;
    z = Optizelle.SQL.create( ...
        [Optizelle.Cone.Semidefinite, ...
         Optizelle.Cone.Quadratic, ...
         Optizelle.Cone.Linear], ...
        [2,2,1]);
    z.data{1} = [
        x(1) x(2);
        x(2) x(3)];
    z.data{2} = [x(2);x(3)];
    z.data{3} = x(3)-1;
end

% z=h'(x)dx
function z=MyIneq_p(x,dx)
    global Optizelle;
    z = Optizelle.SQL.create( ...
        [Optizelle.Cone.Semidefinite, ...
         Optizelle.Cone.Quadratic, ...
         Optizelle.Cone.Linear], ...
        [2,2,1]);
    z.data{1} = [
        dx(1) dx(2);
        dx(2) dx(3)];
    z.data{2} = [dx(2);dx(3)];
    z.data{3} = dx(3);
end

% x_hat=h'(x)*dz
function x_hat=MyIneq_ps(x,dz)
    % Remember, the input to this function may not be symmetric, so
    % compute accordingly 
    x_hat = ...
        [dz.data{1}(1,1); ...
         dz.data{1}(1,2) + dz.data{1}(2,1); ...
         dz.data{1}(2,2)] + ...
        [0;dz.data{2}] + ...
        [0;0;dz.data{3}];
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess for the primal
    x = [5.5; 3.3; 2.2];

    % Generate an initial guess for the dual
    %---SQLSpec0---
    types = ...
        [Optizelle.Cone.Semidefinite, ...
         Optizelle.Cone.Quadratic, ...
         Optizelle.Cone.Linear];
    sizes = [2,2,1];
    %---SQLSpec1---
    %---SQLVector0---
    z = Optizelle.SQL.create(types,sizes);
    %---SQLVector1---

    % Create an optimization state
    state=Optizelle.InequalityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.SQL,x,z);

    % Read the parameters from file
    state=Optizelle.json.InequalityConstrained.read( ...
        Optizelle.Rm,Optizelle.SQL,fname,state);
    
    % Create a bundle of functions
    fns=Optizelle.InequalityConstrained.Functions.t;
    fns.f=MyObj();
    fns.h=MyIneq();

    % Solve the optimization problem
    state=Optizelle.InequalityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging.stdout,fns,state);
    
    % Write an intermediate restart file 
    Optizelle.json.InequalityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.SQL,'restart.json',state);

    % Read in the restart file
    state = Optizelle.json.InequalityConstrained.read_restart( ...
        Optizelle.Rm,Optizelle.SQL,'restart.json',state.x,state.z);

    % Change the maximum number of iterations to something larger and reset
    % the convergence flag
    state.iter_max = 500;
    state.opt_stop = Optizelle.OptimizationStop.NotConverged;

    % Finish solving the optimization problem
    state=Optizelle.InequalityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging.stdout,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e)\n', ...
        state.x(1),state.x(2),state.x(3));

    % Write out the final answer to file
    Optizelle.json.InequalityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.SQL,'solution.json',state);
end
