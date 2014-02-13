% This demonstrates how to solve a simple parameter estimation problem
%
% min .5 || x2 - d ||^2 + .5 beta || x1 || ^2 s.t. (sum_i A_i x1_i)x_2 = b
%
% where A_i, d, and b are randomly generated.  Note, since they're randomly
% generated, it's possible to create a degenerate problem that causes
% problems for our linear solvers.

function parest(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('parest <parameters>');
    end

    % Execute the optimization
    main(fname);
end


% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    Optizelle = setupOptizelle();

    % Set the size of the problem
    m = 2; % Size of x1
    n = 3; % Size of x2

    % Create the data
    d=randn(n,1);

    % Create the rhs
    b=randn(n,1);

    % Create the model
    clear A;
    for i=1:m
        A{i}=randn(n);
    end

    % ff is the expanded objecitve function with two pieces x{1} and x{2}
    beta = 1e-2; % Regularization parameter
    ff.eval=@(x) .5*norm(x{2}-d,2)^2 + .5*beta*norm(x{1},2)^2;
    ff.grad_1=@(x) beta*x{1};
    ff.grad_2=@(x) x{2}-d;
    ff.hessvec_11=@(x,dx) beta*dx;
    ff.hessvec_12=@(x,dx) zeros(m,1);
    ff.hessvec_21=@(x,dx) zeros(n,1);
    ff.hessvec_22=@(x,dx) dx;

    % gg is the expanded constraint where we split up the variable into two
    % pieces x{1} and x{2}.
    gg.eval=@(x) g_eval(A,b,x); 
    gg.p_1=@(x,dx1) gp_x1(A,x,dx1); 
    gg.p_2=@(x,dx2) gp_x2(A,x,dx2); % Need inverse
    gg.p_2_inv=@(x,dx2) gp_x2_inv(A,x,dx2); 
    gg.ps_1=@(x,dy) gps_x1(A,x,dy); 
    gg.ps_2=@(x,dy) gps_x2(A,x,dy); % Need inverse
    gg.ps_2_inv=@(x,dy) gps_x2_inv(A,x,dy); 
    gg.pps_11=@(x,dx,dy) zeros(m,1);
    gg.pps_21=@(x,dx,dy) gps_x2(A,{dx,zeros(n,1)},dy); 
    gg.pps_12=@(x,dx,dy) gps_x1(A,{zeros(m,1),dx},dy); 
    gg.pps_22=@(x,dx,dy) zeros(n,1);

    % Create the equality constrained problem
    [RmxRm f g]=genEqualityConstrained(Optizelle.Rm,Optizelle.Rm,ff,gg);

    % Set a starting guess
    x={randn(m,1),randn(n,1)};
    y=zeros(n,1);

    % Create an optimization state
    state = Optizelle.EqualityConstrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,x,y);

    % Read the parameters from file
    state = Optizelle.json.EqualityConstrained.read( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging,fname,state);

    % Create a bundle of functions
    fns=Optizelle.EqualityConstrained.Functions.t;
    fns.f=f;
    fns.g=g;

    % Solve the full-space problem 
    fprintf('\n------------Solving the full-space problem------------\n');
    state = Optizelle.EqualityConstrained.Algorithms.getMin( ...
        RmxRm,Optizelle.Rm,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));

    % Form the solution operator.  The user must specify this
    phi.eval=@(x)phi_eval(A,b,x);

    % Create the unconstrained problem
    f = genUnconstrained(Optizelle.Rm,Optizelle.Rm,ff,gg,phi);

    % Set a starting guess
    x=randn(m,1);
    
    % Create an unconstrained state 
    state=Optizelle.Unconstrained.State.t(Optizelle.Rm,Optizelle.Messaging,x);

    % Read the parameters from file
    state=Optizelle.json.Unconstrained.read(Optizelle.Rm,Optizelle.Messaging,...
        fname,state);

    % Create a bundle of functions
    fns=Optizelle.Unconstrained.Functions.t;
    fns.f=f;

    % Solve the optimization problem
    state = Optizelle.Unconstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Messaging,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.StoppingCondition.to_string(state.opt_stop));
    end

    % Create the functions for the constraint
    function z = g_eval(A,b,x)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x{1}(i)*A{i};
        end

        % Now, find (sum A_i x1_i)x2 - f
        z=B*x{2}-b;
    end

    function z = gp_x1(A,x,dx1)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i dx1_i
        B = zeros(n);
        for i=1:m
            B = B + dx1(i)*A{i};
        end

        % Now, find (sum A_i dx1_i)y
        z=B*x{2};
    end

    function z = gp_x2(A,x,dx2)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x{1}(i)*A{i};
        end

        % Now, find (sum A_i x1_i)dx2
        z=B*dx2;
    end

    function z = gp_x2_inv(A,x,dx2)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x{1}(i)*A{i};
        end

        % Now, solve (sum A_i x1_i) z = dx2
        z=B\dx2;
    end

    function z = gps_x1(A,x,dy)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);

        % Form the adjoint
        z=zeros(m,1);
        for i=1:m
            z(i)=x{2}'*A{i}'*dy;
        end
    end

    function z = gps_x2(A,x,dy)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x{1}(i)*A{i}';
        end

        % Now, find (sum A_i' x1_i)dy
        z=B*dy;
    end

    function z = gps_x2_inv(A,x,dy)
        % Get the sizes
        m = size(x{1},1);
        n = size(x{2},1);
        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x{1}(i)*A{i}';
        end

        % Now, solve (sum A_i' x1_i) z = dy
        z=B\dy;
    end

    % Create the functions for the solution operator
    function z = phi_eval(A,b,x)
        % Cache the result if possible
        persistent zz;
        persistent xx;

        % Get the sizes
        m = size(x,1);
        n = length(b);

        % Check if we've already evaluated the function.
        if length(xx)==m
            if norm(x-xx)/(1e-16+norm(x)) < 1e-15
                z=zz;
                return;
            end
        end

        
        % First, form sum A_i x1_i
        B = zeros(n);
        for i=1:m
            B = B + x(i)*A{i};
        end

        % Now, solve (sum A_i x_i) z = b 
        z=B\b;

        % Cache the result
        xx=x;
        zz=z;
    end
