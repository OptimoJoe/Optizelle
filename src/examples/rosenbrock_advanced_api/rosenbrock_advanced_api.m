% In this example, we duplicate the Rosenbrock example while demonstrating
% some of the more advanced API features such as custom vector spaces,
% messaging objects, and restarts.

function rosenbrock_advanced_api(pname,rname)
    % Read in the name for the input file
    if ~(nargin==1 || nargin==2)
        error(sprintf('%s\n%s', ...
            'rosenbrock_advanced_api(parameters)\n', ...
            'rosenbrock_advanced_api(parameters,restart)'));
    end

    % Execute the optimization
    if nargin==1
        main(pname);
    else
        main(pname,rname);
    end

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
    
    % Jordan product, z <- x o y.
    self.prod = @(x,y)x.*y;

    % Identity element, x <- e such that x o e = x.
    self.id = @(x)ones(size(x));

    % Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
    self.linv = @(x,y)y./x;
    
    % Barrier function, barr <- barr(x) where x o grad barr(x) = e.
    self.barr = @(x)sum(log(x));
    
    % Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
    % where y > 0. 
    self.srch = @(x,y) feval(@(z)min([min(z(find(z>0)));inf]),-y./x);
    
    % Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
    % operator.
    self.symm = @(x)x;
end

% Squares its input
function z = sq(x)
    z=x*x;
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

% Define a perfect preconditioner for the Hessian
function self = RosenHInv()
    self.eval = @(state,dx) eval(state,dx);

    function result = eval(state,dx)
        x = state.x;
        one_over_det=1./(80000.*sq(x(1))-80000.*x(2)+400.);
        result = [
            one_over_det*(200.*dx(1)+400.*x(1)*dx(2));
            one_over_det*...
                (400.*x(1)*dx(1)+(1200.*x(1)*x(1)-400.*x(2)+2.)*dx(2))];
    end
end

% Define a custom messaging object
function msg = MyMessaging()
    msg = struct( ...
        'print',@(x)fprintf('PRINT:  %s\n',x), ...
        'error',@(x)error(sprintf('ERROR:  %s',x)));
end

% Define serialization routines for MyVS
function MySerialization()
    Optizelle.json.Serialization.serialize( ...
        'register',@(x)strrep(mat2str(x'),' ',', '),@isvector);
    Optizelle.json.Serialization.deserialize( ...
        'register',@(x,x_json)str2num(x_json)',@isvector);
end
    
% Define a state manipulator that writes out the optimization state at
% each iteration.
function smanip=MyRestartManipulator()
    smanip=struct('eval',@(fns,state,loc)MyRestartManipulator_(fns,state,loc));
end
function state=MyRestartManipulator_(fns,state,loc)
    global Optizelle;

    % At the end of the optimization iteration, write the restart file
    if(loc == Optizelle.OptimizationLocation.EndOfOptimizationIteration)
        % Create a reasonable file name
        ss = sprintf('rosenbrock_advanced_api_%04d.json',state.iter);
            
        % Write the restart file
        Optizelle.json.Unconstrained.write_restart( ...
           MyVS(),MyMessaging(),ss,state);
    end
end

% Actually runs the program
function main(pname,rname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess for Rosenbrock
    x = [-1.2;1.];

    % Create an unconstrained state based on this vector
    state=Optizelle.Unconstrained.State.t(MyVS(),MyMessaging(),x);
    
    % If we have a restart file, read in the parameters 
    if(nargin==2)
        state = Optizelle.json.Unconstrained.read_restart( ...
            MyVS(),MyMessaging(),rname,x,state);
    end
    
    % Read the parameters from file
    state=Optizelle.json.Unconstrained.read( ...
        MyVS(),MyMessaging(),pname,state);

    % Create the bundle of functions 
    fns=Optizelle.Unconstrained.Functions.t;
    fns.f=Rosenbrock();
    fns.PH=RosenHInv();

    % Solve the optimization problem
    state=Optizelle.Unconstrained.Algorithms.getMin( ...
        MyVS(),MyMessaging(),MyRestartManipulator(),fns,state);
    
    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x(1),state.x(2));
    
    % Write out the final answer to file
    Optizelle.json.Unconstrained.write_restart( ...
        MyVS(),MyMessaging(),'solution.json',state);
end
