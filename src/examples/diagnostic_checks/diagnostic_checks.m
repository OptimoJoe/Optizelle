% This example demonstrates how to run a series of diagnostic tests
% on functions and then immediately exit.
function diagnostic_checks()
    % Execute the optimization
    main();
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Cubes its input
function z = cub(x)
    z=x*x*x;
end

% Quads its input
function z = quad(x)
    z=x*x*x*x;
end

% Quints its input
function z = quint(x)
    z=x*x*x*x*x;
end

% Define the Rosenbrock function where
% 
% f(x,y)=(1-x)^2+100(y-x^2)^2
%
function self = Rosenbrock()

    % Evaluation of the Rosenbrock function
    self.eval = @(x) sq(1.-x(1))+100.*sq(x(2)-sq(x(1)));

    % Gradient
    self.grad = @(x) [
        -400.*x(1)*(x(2)-sq(x(1)))-2.*(1.-x(1));
        200.*(x(2)-sq(x(1)))];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
        (1200.*sq(x(1))-400.*x(2)+2)*dx(1)-400.*x(1)*dx(2);
        -400.*x(1)*dx(1)+200.*dx(2)];
end

% Define some utility function where
%
% g(x)= [ cos(x1) sin(x2)   ]
%       [ 3 x1^2 x2 + x2 ^3 ]
%       [ log(x1) + 3 x2 ^5 ]
%
function self=Utility()

    % y=g(x) 
    self.eval = @(x) [
        cos(x(1))*sin(x(2));
        3.*sq(x(1))*x(2)+cub(x(2));
        log(x(1))+3.*quint(x(2))];

    % y=g'(x)dx
    self.p = @(x,dx) [
        ( -sin(x(1))*sin(x(2))*dx(1) ...
          +cos(x(1))*cos(x(2))*dx(2));
        ( 6.*x(1)*x(2)*dx(1) ...
          +(3.*sq(x(1))+3.*sq(x(2)))*dx(2));
        ( 1./x(1)*dx(1) ...
          +15.*quad(x(2))*dx(2))];

    % z=g'(x)*dy
    self.ps = @(x,dy) [
        ( -sin(x(1))*sin(x(2))*dy(1) ...
          +6.*x(1)*x(2)*dy(2) ...
          +1./x(1)*dy(3));
        ( cos(x(1))*cos(x(2))*dy(1) ...
          +(3.*sq(x(1))+3.*sq(x(2)))*dy(2) ...
          +15.*quad(x(2))*dy(3))];

    % z=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) [
        (  (-cos(x(1))*dx(1)*sin(x(2))-sin(x(1))*cos(x(2))*dx(2))*dy(1) ...
          +(6.*dx(1)*x(2) + 6.*x(1)*dx(2))*dy(2) ...
          +(-1./sq(x(1))*dx(1))*dy(3));
        (  (-sin(x(1))*dx(1)*cos(x(2))-cos(x(1))*sin(x(2))*dx(2))*dy(1) ...
          +(6.*x(1)*dx(1)+6.*x(2)*dx(2))*dy(2) ...
          +(60.*cub(x(2))*dx(2))*dy(3)) ];
end

% Actually runs the program
function main(fname)
    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle ();

    % Allocate memory for an initial guess and equality multiplier 
    x = [1.2;2.3];
    z = zeros(3,1);

    % Create an optimization state
    state=Optizelle.EqualityConstrained.State.t(Optizelle.Rm,Optizelle.Rm,x,z);

    %Modify the state so that we just run our diagnostics and exit
    state.dscheme = Optizelle.DiagnosticScheme.DiagnosticsOnly;
    state.f_diag = Optizelle.FunctionDiagnostics.SecondOrder;
    state.x_diag = Optizelle.VectorSpaceDiagnostics.Basic;
    state.h_diag = Optizelle.FunctionDiagnostics.SecondOrder;
    state.z_diag = Optizelle.VectorSpaceDiagnostics.EuclideanJordan;
    state.L_diag = Optizelle.FunctionDiagnostics.SecondOrder;

    % Create a bundle of functions
    fns=Optizelle.EqualityConstrained.Functions.t;
    fns.f=Rosenbrock();
    fns.g=Utility();

    % Even though this looks like we're solving an optimization problem,
    % we're actually just going to run our diagnostics and then exit.
    Optizelle.EqualityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);
end
