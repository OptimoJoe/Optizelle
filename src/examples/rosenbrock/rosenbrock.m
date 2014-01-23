% In this example, we setup and minimize the Rosenbrock function.
function rosenbrock(fname)
    % Read in the name for the input file
    if nargin ~=2
        error('rosenbrock <parameters>');
    end

    % Execute the optimization
    main(fname);
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Define the Rosenbrock function where
% 
% f(x,y)=(1-x)^2+100(y-x^2)^2
%
function self = Rosenbrock(self)
    
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

% Define a perfect preconditioner for the Hessian
function self = RosenHInv(self)
    self.eval = @(state,dx) eval(state,dx);

    function result = eval(state,dx)
        x = state.x;
        double one_over_det=1./(400000.*x(1)*x(1)-80000.*x(2)+400.);
        result = [
            one_over_det*(200.*dx(1)+400.*x(1)*dx(2));
            one_over_det*...
                (400.*x(1)*dx(1)+(1200.*x(1)*x(1)-400.*x(2)+2.)*dx(2))];
    end
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    Optizelle = setupOptizelle ();

    % Generate an initial guess for Rosenbrock
    x = [-1.2;1.];

    % Create an unconstrained state based on this vector
    state=Optizelle.Unconstrained.State.t(Optizelle.Rm,Optizelle.Messaging,x);

    % Read the parameters from file
    Optizelle.json.Unconstrained.read(Optizelle.Rm,Optizelle.Messaging, ...
        fname,state);

    % Create the bundle of functions 
    fns.Optizelle.Unconstrained.Functions.t;
    fns.f=Rosenrock(Optizelle.ScalarValuedFunction);
    fns.PH=RosenHInv(Optizelle.Operator);
    
    % Solve the optimization problem
    Optizelle.Unconstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));

    % Write out the final answer to file
    Optizelle.json.Unconstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Messaging,'solution.json',state);
end
