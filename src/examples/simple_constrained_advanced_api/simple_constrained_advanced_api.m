% Optimize a simple optimization problem with an optimal solution
% of (1/3,1/3)

function simple_constrained_advanced_api(pname,rname)
    % Read in the name for the input file
    if ~(nargin==1 || nargin==2)
        error(sprintf('%s\n%s', ...
            'simple_constrained_advanced_api(parameters)\n', ...
            'simple_constrained_advanced_api(parameters,restart)'));
    end

    % Execute the optimization
    if nargin==1
        main(pname);
    else
        main(pname,rname);
    end
end

% Convert a vector to structure
function y = tostruct(x)
    y = struct('data',x);
end

% Defines the vector space used for optimization.
function self = MyVS()

    % Memory allocation and size setting
    self.init = @(x) x;

    % <- x (Shallow.  No memory allocation.)
    self.copy = @(x) x;

    % <- alpha * x
    self.scal = @(alpha,x) tostruct(alpha*x.data);

    % <- 0
    self.zero = @(x) tostruct(zeros(size(x.data)));

    % <- alpha * x + y
    self.axpy = @(alpha,x,y) tostruct(alpha * x.data + y.data);

    %<- <x,y>
    self.innr = @(x,y)x.data'*y.data;

    % <- random
    self.rand = @(x)tostruct(randn(size(x.data)));
    
    % Jordan product, z <- x o y.
    self.prod = @(x,y)tostruct(x.data .* y.data);

    % Identity element, x <- e such that x o e = x.
    self.id = @(x)tostruct(ones(size(x.data)));

    % Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
    self.linv = @(x,y)tostruct(y.data ./ x.data);
    
    % Barrier function, barr <- barr(x) where x o grad barr(x) = e.
    self.barr = @(x)sum(log(x.data));
    
    % Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
    % where y > 0. 
    self.srch = @(x,y) feval(@(z)min([min(z(find(z>0)));inf]),-y.data ./x.data);
    
    % Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
    % operator.
    self.symm = @(x)x;
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Define a simple objective where 
% 
% f(x,y)=(x+1)^2+(y+1)^2
%
function self = MyObj()

    % Evaluation 
    self.eval = @(x) feval(@(x)sq(x(1)+1.)+sq(x(2)+1.),x.data);

    % Gradient
    self.grad = @(x) tostruct(feval(@(x)[
        2.*x(1)+2.;
        2.*x(2)+2.],x.data));

    % Hessian-vector product
    self.hessvec = @(x,dx) tostruct(feval(@(x,dx)[
        2.*dx(1);
        2.*dx(2)],x.data,dx.data));
end

% Define a simple equality
%
% g(x,y)= [ x + 2y = 1 ] 
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) tostruct(feval(@(x)[x(1)+2.*x(2)-1.],x.data));

    % y=g'(x)dx
    self.p = @(x,dx) tostruct(feval(@(x,dx)[dx(1)+2.*dx(2)],x.data,dx.data));

    % z=g'(x)*dy
    self.ps = @(x,dy) tostruct(feval(@(x,dy)[
        dy(1);
        2.*dy(1)],x.data,dy.data));

    % z=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) tostruct(zeros(2,1)); 
end

% Define simple inequalities 
%
% h(x,y)= [ 2x + y >= 1 ] 
%
function self = MyIneq()

    % y=h(x) 
    self.eval = @(x) tostruct(feval(@(x)[
        2.*x(1)+x(2)-1],x.data));

    % y=h'(x)dx
    self.p = @(x,dx) tostruct(feval(@(x,dx)[
        2.*dx(1)+dx(2)],x.data,dx.data));

    % z=h'(x)*dy
    self.ps = @(x,dy) tostruct(feval(@(x,dy)[
        2.*dy(1)
        dy(1)],x.data,dy.data));

    % z=(h''(x)dx)*dy
    self.pps = @(x,dx,dy) tostruct([ 0. ]); 
end

% Define the serialize routine for MyVS
function x_json=serialize_MyVS(x,name,iter)
    % Create the filename where we put our vector
    fname=sprintf('./restart/%s.%04d.txt',name,iter);

    % Actually write the vector there
    dlmwrite(fname,x.data);

    % Use this filename as the json string 
    x_json = sprintf('\"%s\"',fname);
end

% Define the deserialize routine for MyVS
function x=deserialize_MyVS(x_,x_json)
    % Filter out the quotes and newlines from the string
    x_json = strrep(x_json,'"','');
    x_json = strrep(x_json,sprintf('\n'),'');

    % Read the data into x
    x=tostruct(dlmread(x_json));
end

% Define serialization routines for MyVS
function MySerialization()
    global Optizelle;
    Optizelle.json.Serialization.serialize( ...
        'register', ...
        @(x,name,iter)serialize_MyVS(x,name,iter), ...
        @(x)isstruct(x) && isfield(x,'data') && isvector(x.data));
    Optizelle.json.Serialization.deserialize( ...
        'register', ...
        @(x,x_json)deserialize_MyVS(x,x_json), ...
        @(x)isstruct(x) && isfield(x,'data') && isvector(x.data));
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
        ss = sprintf('simple_constrained_advanced_api_%04d.json',state.iter);
            
        % Write the restart file
        Optizelle.json.Constrained.write_restart( ...
           MyVS(),MyVS(),MyVS(),Optizelle.Messaging,ss,state);
    end
end

% Actually runs the program
function main(pname,rname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Register the serialization routines
    MySerialization();

    % Generate an initial guess 
    x = tostruct([2.1;1.1]);

    % Allocate memory for the equality multiplier 
    y = tostruct([0.]);

    % Allocate memory for the inequality multiplier 
    z = tostruct([0.]);

    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        MyVS(),MyVS(),MyVS(),Optizelle.Messaging,x,y,z);
    
    % If we have a restart file, read in the parameters 
    if(nargin==2)
        state = Optizelle.json.Constrained.read_restart( ...
            MyVS(),MyVS(),MyVS(),Optizelle.Messaging,rname,x,y,z);
    end

    % Read the parameters from file
    state = Optizelle.json.Constrained.read( ...
        MyVS(),MyVS(),MyVS(),Optizelle.Messaging,pname,state);

    % Create a bundle of functions
    fns = Optizelle.Constrained.Functions.t;
    fns.f = MyObj();
    fns.g = MyEq();
    fns.h = MyIneq();

    % Solve the optimization problem
    state = Optizelle.Constrained.Algorithms.getMin( ...
        MyVS(),MyVS(),MyVS(),Optizelle.Messaging,fns,state, ...
        MyRestartManipulator());

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e)\n',state.x.data(1),state.x.data(2));

    % Write out the final answer to file
    Optizelle.json.Constrained.write_restart( ...
        MyVS(),MyVS(),MyVS(), ...
        Optizelle.Messaging,'solution.json',state);
end
