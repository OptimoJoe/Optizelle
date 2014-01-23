% This example minimizes the Rosenbrock function from scratch.  Meaning,
% it runs through a complete example from defining a valid vector space,
% to setting parameters, to solving the problem.  For reference, the optimal
% solution to the Rosenbrock function is (1,1).
function rosenbrock_custom_VS(fname)
    % Read in the name for the input file
    if nargin ~=2
        error('rosenbrock_custom_VS <parameters>');
    end

    % Execute the optimization
    main(fname);
end


% Defines the vector space used for optimization.
function self = MyVS()
    % Memory allocation and size setting
    self.init = @(x) x;

    % <- x (Shallow.  No memory allocation.)
    self.copy = @(x) x;

    % <- alpha * x
    self.scal = @(alpha,x) alpha*x;

    % <- 0
    self.zero = @(x) zeros(size(x));

    % <- alpha * x + y
    self.axpy = @(alpha,x,y) alpha * x + y;

    %<- <x,y>
    self.innr = @(x,y)x'*y;

    % <- random
    self.rand = @(x)randn(size(x));
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

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    Optizelle = setupOptizelle ();

    % Generate an initial guess for Rosenbrock
    x = [-1.2;1.];

    % Create an unconstrained state based on this vector
    state=Optizelle.Unconstrained.State.t(MyVS(),Optizelle.Messaging(),x);

    % Setup some algorithmic parameters

    if 1,
        % Trust-Region Newton's method
        state.H_type = Optizelle.Operators.UserDefined;
        state.iter_max = 50;
        state.eps_krylov = 1e-10;
    end

    if 0,
        % BFGS
        state.algorithm_class = Optizelle.AlgorithmClass.LineSearch;
        state.dir = Optizelle.LineSearchDirection.BFGS;
        state.stored_history = 10;
        state.iter_max = 100;
    end

    if 0,
        % Newton-CG 
        state.algorithm_class = Optizelle.AlgorithmClass.LineSearch;
        state.dir = Optizelle.LineSearchDirection.NewtonCG;
        state.H_type = Optizelle.Operators.UserDefined;
        state.eps_krylov = 1e-2;
        state.iter_max = 50;
    end

    % Create the bundle of functions 
    fns=Optizelle.Unconstrained.Functions.t;
    fns.f=Rosenbrock()

    % Solve the optimization problem
    Optizelle.Unconstrained.Algorithms.getMin(MyVS(),Optizelle.Messaging, ...
        fns,state)

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));
end
