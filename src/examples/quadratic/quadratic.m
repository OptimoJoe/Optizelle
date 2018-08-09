% This example minimizes a simple quadratic function.  For reference, the
% optimal solution to this function is (1,1,1).
function quadratic(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('quadratic <parameters>')
    end

    % Execute the optimization
    main(fname);
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Define the quadratic function
%
% f(x,y,z)=(x-1)^2+(2y-2)^2+(3z-3)^2
%
function self = Quad()

    % Evaluation of the quadratic function
    self.eval = @(x) sq(x(1)-1.)+sq(2*x(2)-2.)+sq(3*x(3)-3.);

    % Gradient
    self.grad = @(x) [
        2*x(1)-2;
        8*x(2)-8;
        18*x(3)-18];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
    	2*dx(1);
        8*dx(2);
        18*dx(3)];
end

% Define an almost perfect preconditioner for the Hessian
function self = QuadHInv()
    self.eval = @(state,dx) [
        dx(1)/2.;
        dx(2)/8.;
        dx(3)];
end


% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess
    x = [-1.2;1.1;2.];

    % Create an unconstrained state based on this vector
    state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x);

    % Read the parameters from file
    state=Optizelle.json.Unconstrained.read(Optizelle.Rm,fname,state);

    % Create the bundle of functions
    fns=Optizelle.Unconstrained.Functions.t;
    fns.f=Quad();
    fns.PH=QuadHInv();

    % Solve the optimization problem
    state=Optizelle.Unconstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e)\n', ...
        state.x(1),state.x(2),state.x(3));

    % Write out the final answer to file
    Optizelle.json.Unconstrained.write_restart( ...
        Optizelle.Rm,'solution.json',state);
end
