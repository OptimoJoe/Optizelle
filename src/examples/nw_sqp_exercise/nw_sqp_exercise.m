% Exercise 18.2 in Numerical Optimization by Nocedal and Wright.  This has
% an optimal solution of x = (-1.71,1.59,1.82.-0.763,-0.763).
function nw_sqp_exercise(fname)
    % Read in the name for the input file
    if nargin ~=1
        error('nw_sqp_exercise <parameters>')
    end

    % Execute the optimization
    main(fname);
end

% Squares its input
function z = sq(x)
    z=x*x;
end

% Power function
function z = pow(x,n)
   z=x^n; 
end

% Indexing for vectors
function result = itok(i)
    result = i;
end

% Indexing for packed storage
function result = ijtokp(i,j)
    if i>j
        tmp=i;
        i=j;
        j=tmp;
    end
    result = i+j*(j-1)/2;
end
    
% Indexing function for dense matrices 
function result = ijtok(i,j,m)
    result = i+(j-1)*m;
end

% Indexing function for dense tensors 
function result = ijktol(i,j,k,m,n)
    result = i+(j-1)*m+(k-1)*m*n;
end

%
% f(x) = exp(x1 x2 x3 x4 x5) - (1/2) (x1^3 + x2^3 + 1)^2
%
function self = MyObj()

    % Evaluation 
    self.eval = @(x) ...
        (exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))) ...
            - sq(pow(x(itok(1)),3)+pow(x(itok(2)),3)+1.)/2.);

    % Gradient
    self.grad= @(x) [ ...
        (x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))) ...
            - 3.*sq(x(itok(1))) ...
            *(pow(x(itok(1)),3) + pow(x(itok(2)),3)+1.));
        (x(itok(1))*x(itok(3))*x(itok(4))*x(itok(5)) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))) ...
            - 3.*sq(x(itok(2))) ...
            * (pow(x(itok(1)),3) + pow(x(itok(2)),3) + 1.));
        (x(itok(1))*x(itok(2))*x(itok(4))*x(itok(5)) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))));
        (x(itok(1))*x(itok(2))*x(itok(3))*x(itok(5)) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))));
        (x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4)) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))))];

    % Hessian-vector product
    self.hessvec = @(x,dx) hessvec(x,dx);
end
function H_dx = hessvec(x,dx)
    % Allocate memory for the dense Hessian in packed storage
    H=zeros(15,1);

    % Compute the dense Hessian
    H(ijtokp(1,1)) = ( ...
        sq(x(itok(2)))*sq(x(itok(3)))*sq(x(itok(4)))*sq(x(itok(5))) ...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            - 6.*x(itok(1))*(pow(x(itok(1)),3) ...
            + pow(x(itok(2)),3)+1.) ...
            - 9.*pow(x(itok(1)),4));
    H(ijtokp(1,2)) = ( ...
       x(itok(1))*x(itok(2))*sq(x(itok(3)))*sq(x(itok(4)))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(3))*x(itok(4))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5)))...
            - 9.*sq(x(itok(1)))*sq(x(itok(2))));
    H(ijtokp(1,3)) = ( ...
       x(itok(1))*sq(x(itok(2)))*x(itok(3))*sq(x(itok(4)))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(2))*x(itok(4))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(1,4)) = ( ...
       x(itok(1))*sq(x(itok(2)))*sq(x(itok(3)))*x(itok(4))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(2))*x(itok(3))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(1,5)) = ( ...
       x(itok(1))*sq(x(itok(2)))*sq(x(itok(3)))*sq(x(itok(4)))*x(itok(5))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(2))*x(itok(3))*x(itok(4))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(2,2)) = ( ...
        sq(x(itok(1)))*sq(x(itok(3)))*sq(x(itok(4)))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            - 6.*x(itok(2)) ...
            * (pow(x(itok(1)),3)+pow(x(itok(2)),3)+1.) ...
            - 9.*pow(x(itok(2)),4));
    H(ijtokp(2,3)) = ( ...
       sq(x(itok(1)))*x(itok(2))*x(itok(3))*sq(x(itok(4)))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(4))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(2,4))= ( ...
       sq(x(itok(1)))*x(itok(2))*sq(x(itok(3)))*x(itok(4))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(3))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(2,5))= ( ...
       sq(x(itok(1)))*x(itok(2))*sq(x(itok(3)))*sq(x(itok(4)))*x(itok(5))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(3))*x(itok(4))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(3,3)) = ( ...
        sq(x(itok(1)))*sq(x(itok(2)))*sq(x(itok(4)))*sq(x(itok(5))) ...
        * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(3,4)) = ( ...
       sq(x(itok(1)))*sq(x(itok(2)))*x(itok(3))*x(itok(4))*sq(x(itok(5)))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(2))*x(itok(5))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(3,5)) = ( ...
       sq(x(itok(1)))*sq(x(itok(2)))*x(itok(3))*sq(x(itok(4)))*x(itok(5))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(2))*x(itok(4))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(4,4)) = ( ...
        sq(x(itok(1)))*sq(x(itok(2)))*sq(x(itok(3)))*sq(x(itok(5))) ...
        * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))); 
    H(ijtokp(4,5)) = ( ...
       sq(x(itok(1)))*sq(x(itok(2)))*sq(x(itok(3)))*x(itok(4))*x(itok(5))...
            * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5)))...
            + x(itok(1))*x(itok(2))*x(itok(3))*exp(x(itok(1))*x(itok(2))...
            * x(itok(3))*x(itok(4))*x(itok(5))));
    H(ijtokp(5,5)) = ( ...
        sq(x(itok(1)))*sq(x(itok(2)))*sq(x(itok(3)))*sq(x(itok(4))) ...
        * exp(x(itok(1))*x(itok(2))*x(itok(3))*x(itok(4))*x(itok(5))));

    % Compute the Hessian-vector product
    H_dx = zeros(5,1);
    for i = 1:5 
        for j = 1:5 
            H_dx(i) = H_dx(i) + H(ijtokp(i,j))*dx(j);
        end
    end
end

%
% g(x)= [ x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 ]
%       [ x2 x3 - 5 x4 x5                       ]
%       [ x1^3 + x2^3 + 1                       ]
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) [
        (sq(x(itok(1))) + sq(x(itok(2))) + sq(x(itok(3))) ...
            + sq(x(itok(4))) + sq(x(itok(5))) - 10.);
        (x(itok(2))*x(itok(3)) - 5.*x(itok(4))*x(itok(5)));
        (pow(x(itok(1)),3) + pow(x(itok(2)),3) + 1.)];

    % y=g'(x)dx
    self.p = @(x,dx) reshape(generateJac(x),3,5)*dx; 

    % xhat=g'(x)*dy
    self.ps = @(x,dy) reshape(generateJac(x),3,5)'*dy; 

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) pps(x,dx,dy); 
end
% Generate a dense version of the Jacobian
function jac = generateJac(x)
    jac = zeros(15,1);

    jac(ijtok(1,1,3)) = 2.*x(itok(1));
    jac(ijtok(1,2,3)) = 2.*x(itok(2));
    jac(ijtok(1,3,3)) = 2.*x(itok(3));
    jac(ijtok(1,4,3)) = 2.*x(itok(4));
    jac(ijtok(1,5,3)) = 2.*x(itok(5));
    
    jac(ijtok(2,2,3)) = x(itok(3));
    jac(ijtok(2,3,3)) = x(itok(2));
    jac(ijtok(2,4,3)) = -5.*x(itok(5));
    jac(ijtok(2,5,3)) = -5.*x(itok(4));
    
    jac(ijtok(3,1,3)) = 3.*sq(x(itok(1)));
    jac(ijtok(3,2,3)) = 3.*sq(x(itok(2)));
end
function z = pps(x,dx,dy)
    % Generate a dense tensor that holds the second derivative adjoint
    D = zeros(75,1);
    D(ijktol(1,1,1,3,5)) = 2.;
    D(ijktol(1,2,2,3,5)) = 2.;
    D(ijktol(1,3,3,3,5)) = 2.;
    D(ijktol(1,4,4,3,5)) = 2.;
    D(ijktol(1,5,5,3,5)) = 2.;
    
    D(ijktol(2,2,3,3,5)) = 1.;
    D(ijktol(2,3,2,3,5)) = 1.;
    D(ijktol(2,4,5,3,5)) = -5.;
    D(ijktol(2,5,4,3,5)) = -5.;

    D(ijktol(3,1,1,3,5)) = 6.*x(itok(1));
    D(ijktol(3,2,2,3,5)) = 6.*x(itok(2));

    % Compute the action of this operator on our directions
    z=zeros(5,1);
    for i = 1:3,
        for j = 1:5,
            for k = 1:5,
                z(itok(k)) = z(itok(k)) + ...
                    D(ijktol(i,j,k,3,5))*dx(itok(j))*dy(itok(i));
            end
        end
    end
end


% Actually runs the program
function main(fname)
    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess for the primal
    x = [-1.8;1.7;1.9;-0.8;-0.8];

    % Generate an initial guess for the dual
    y = zeros(3,1);

    % Create an optimization state
    state=Optizelle.EqualityConstrained.State.t(Optizelle.Rm,Optizelle.Rm,x,y);

    % Read the parameters from file
    state=Optizelle.json.EqualityConstrained.read(Optizelle.Rm,Optizelle.Rm, ...
        fname,state);

    % Create the bundle of functions 
    fns=Optizelle.EqualityConstrained.Functions.t;
    fns.f=MyObj();
    fns.g=MyEq();

    % Solve the optimization problem
    state = Optizelle.EqualityConstrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is:\n');
    for i=1:5, 
        if i==1,
            fprintf('[ ');
        else
            fprintf('  ');
        end
        fprintf('%13e',state.x(itok(i)));
        if i==5,
            fprintf(' ]\n');
        else
            fprintf(' ;\n');
        end
    end

    % Write out the final answer to file
    Optizelle.json.EqualityConstrained.write_restart( ...
        Optizelle.Rm,Optizelle.Rm,'solution.json',state);
end
