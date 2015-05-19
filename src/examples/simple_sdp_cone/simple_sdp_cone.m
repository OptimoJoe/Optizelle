% Optimize a simple problem with an optimal solution of (0.5,.25)

function simple_sdp_cone(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('simple_sdp_cone <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Define a simple objective where 
% 
% f(x,y)=-x+y
%
function self = MyObj()

    % Evaluation 
    self.eval = @(x) -x(1)+x(2);

    % Gradient
    self.grad = @(x) [
        -1.;
        1.];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
        0;
        0];
end

% Define a simple SDP inequality 
%
% h(x,y) = [ y x ] >= 0
%          [ x 1 ]
%
function self = MyIneq()

    % y=h(x) 
    self.eval = @(x)MyIneq_eval(x);

    % y=h'(x)dx
    self.p = @(x,dx)MyIneq_p(x,dx);

    % z=h'(x)*dy
    self.ps = @(x,dy) [
        dy.data{1}(2,1)+dy.data{1}(1,2);
        dy.data{1}(1,1)];

    % z=(h''(x)dx)*dy
    self.pps = @(x,dx,dy) [
        0;
        0];
end

% y=h(x) 
function y=MyIneq_eval(x)
    global Optizelle;
    y = Optizelle.SQL.create([Optizelle.Cone.Semidefinite],[2]);
    y.data{1} = [
        x(2) x(1);
        x(1) 1.  ];
end

% y=h'(x)dx
function y=MyIneq_p(x,dx)
    global Optizelle;
    y = Optizelle.SQL.create([Optizelle.Cone.Semidefinite],[2]);
    y.data{1} = [
        dx(2) dx(1);
        dx(1) 0.   ];
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess for the primal
    x = [1.2; 3.1];

    % Generate an initial guess for the dual
    z = Optizelle.SQL.create([Optizelle.Cone.Semidefinite],[2]);

    % Create an optimization state
    state=Optizelle.InequalityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging,x,z);

    % Read the parameters from file
    state=Optizelle.json.InequalityConstrained.read( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging,fname,state);
    
    % Create a bundle of functions
    fns=Optizelle.InequalityConstrained.Functions.t;
    fns.f=MyObj();
    fns.h=MyIneq();

    % Solve the optimization problem
    state=Optizelle.InequalityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));

    % Write out the final answer to file
    Optizelle.json.InequalityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.SQL,Optizelle.Messaging,'solution.json',state);
end
