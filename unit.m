%% Set the random seed
randn('seed',0);

% Create some random data
m=10; % size of system response 
n=20; % number of reponses
d=randn(m,n);

% Create a random point for evaluation
y=randn(m,1);

% Check the objective function
for i=1:n
    if norm(f(i,d,y)-(y-d(:,i))) > 1e-12
    	fprintf('Error in %dth data point during evaluation of f\n',i)
	return
    end
end
fprintf('Evaluation of f complete.\n');

% Check the objective function's derivative
eta=randn(m,1);
for i=1:n
    if norm(fp(i,d,y,eta)-eta) > 1e-12
    	fprintf('Error in %dth data point during evaluation of fp\n',i)
	return
    end
end
fprintf('Evaluation of fp complete.\n');

% Check the objective function's derivative adjoint
xi=randn(m,1);
for i=1:n
    if norm(fps(i,d,y,xi)-xi) > 1e-12
    	fprintf('Error in %dth data point during evaluation of fps\n',i)
	return
    end
end
fprintf('Evaluation of fps complete.\n');

% Check the objective function's 2nd derivative adjoint
eta=randn(m,1);
xi=randn(m,1);
for i=1:n
    if norm(fpps(i,d,y,eta,xi)-zeros(m,1)) > 1e-12
    	fprintf('Error in %dth data point during evaluation of fps\n',i)
	return
    end
end
fprintf('Evaluation of fpps complete.\n');

% Create data to test the solution operator 
% Create the linear operators
p=25; % The number of linear operators
A={};
for i=1:p
    A={A{:} randn(m)};
end
B=randn(m);

% Create the data
q=15; % The number of linear system solves 
b={};
for i=1:q
    b={b{:} randn(m,1)};
end

% Create the control
u=randn(1,p);

% Create the complete operator
C=zeros(m);
for i=1:p
    C=C+A{i}*u(i);
end
C=C+B;

% Check the solution operator
for i=1:q
    if norm(h(A,B,b,i,u)-C\b{i}) > 1e-12
    	fprintf('Error in the %dth data point during evaluation of h\n',i)
    end
end
fprintf('Evaluation of h complete.\n');


% Check the derivative of the solution operator
eta=randn(1,p);

% Build the operator -A(eta)
Aeta=zeros(m);
for i=1:p
    Aeta=Aeta-A{i}*eta(i);
end

% Check the derivative
for i=1:q
    if norm(hp(A,B,b,i,u,eta)-(C\(Aeta*h(A,B,b,i,u)))) > 1e-11
    	fprintf('Error in the %dth data point during evaluation of hp\n',i)
    end
end
fprintf('Evaluation of hp complete.\n');

% Check the adjoint of the derivative operator
eta=randn(m,1);

% Do an adjoint solve
xi=C'\eta;

% Check for each rhs 
for i=1:q
    % Find the solution operator
    y=h(A,B,b,i,u);

    % Construct the application of the adjoint on xi
    z=zeros(1,p);
    for j=1:p
    	z(j)=-xi'*A{j}*y;
    end

    % Check the difference
    if norm(hps(A,B,b,i,u,eta)-z) > 1e-12
    	fprintf('Error in the %dth data point during evaluation of hps\n',i)
    end
end
fprintf('Evaluation of hps complete.\n');

% Check the adjoint of the second derivative
eta=randn(1,p);
xi=randn(m,1);

% Do an adjoint solve
nu=C'\xi;

% Check for each rhs
for i=1:q
    % Find Aeta
    Aeta=zeros(m);
    for j=1:p
	Aeta=Aeta+A{j}*eta(j);
    end

    % Find (Aeta*)nu
    Aetat_nu=Aeta'*nu;

    % Find -(h'(u)*)(Aeta*)nu
    tmp1=-hps(A,B,b,i,u,Aetat_nu);

    % Find h'(u)eta
    z=hp(A,B,b,i,u,eta);

    % Find -((A(.)h'(u)eta)*)nu
    tmp2=zeros(1,p);
    for j=1:p
    	tmp2(j)=-nu'*A{j}*z;
    end

    % Take -tmp1-tmp2 to find the solution
    w=tmp1+tmp2;

    % Check the difference
    if norm(hpps(A,B,b,i,u,eta,xi)-w) > 1e-11
    	fprintf('Error in the %dth data point during evaluation of hpps\n',i)
    end
end
fprintf('Evaluation of hpps complete.\n');

% Check the gradient operator

% Create some random data that matches the number of right hand sides
d=randn(m,q);

% Find the gradient
g=getGradient(A,B,b,d,u);

% Find the gradient the long way
gg=zeros(size(u));

% Accumulate each piece of the gradient
for i=1:size(d,2)
    % Find the solution given this control 
    y=h(A,B,b,i,u);

    % Determine how well this matches the data
    ymd=f(i,d,y);

    % Apply f'(y)^* to y-d
    fps_ymd=fps(i,d,y,ymd);

    % apply h'(u)^* to f'(u)^*(y-d)
    gg = gg + hps(A,B,b,i,u,fps_ymd);
end

% Check the difference
if norm(gg-g) > 1e-12
    fprintf('Error when evaluating getGradient\n');
end
fprintf('Evaluation of getGradient complete.\n');

% Check the GN Hessian approximation

% Create a random direction for which we can evaluate
p=randn(size(u));

% Evaluate the GN Hessian approximation in the direction p
q=GaussNewton(A,B,b,d,u,p);

% Find the result the long way
qq=zeros(size(u));

% Accumulate each piece of the GN approximation
for i=1:size(d,2)
    % Find the solution given the control
    y=h(A,B,b,i,u);

    % Find h'(u)p
    hpu_p=hp(A,B,b,i,u,p);

    % Find f'(y)h'(u)p
    fpy_hpu_p=fp(i,d,y,hpu_p);

    % Find f'(y)^*f'(y)h'(u)p
    fpsy_hpy_hpu_p=fps(i,d,y,fpy_hpu_p);

    % Find h'(u)^*f'(y)^*f'(y)h'(u)p
    qq=qq+hps(A,B,b,i,u,fpsy_hpy_hpu_p);
end

% Check the difference
if norm(qq-q) > 1e-12
    fprintf('Error when evaluating GaussNewton\n');
end
fprintf('Evaluation of GaussNewton complete.\n');

% Check the Newton Hessian 

% Create a random direction for which we can evaluate
p=randn(size(u));

% Evaluate the Newton Hessian in the direction p
q=Newton(A,B,b,d,u,p);

% Find the result the long way
qq=zeros(size(u));

% Accumulate each piece of the GN approximation
for i=1:size(d,2)
    % Find the solution given the control
    y=h(A,B,b,i,u);

    % Find h'(u)p
    hpu_p=hp(A,B,b,i,u,p);

    % Find f'(y)h'(u)p
    fpy_hpu_p=fp(i,d,y,hpu_p);

    % Find f'(y)^*f'(y)h'(u)p
    fpsy_hpy_hpu_p=fps(i,d,y,fpy_hpu_p);

    % Find h'(u)^*f'(y)^*f'(y)h'(u)p
    qq=qq+hps(A,B,b,i,u,fpsy_hpy_hpu_p);

    % Find f(y)
    fy=f(i,d,y);

    % Find f'(y)*f(y)
    fpy_fy=fp(i,d,y,fy);

    % Find (h''(u)p)*f'(y)^*f(y)
    hpps_up_fpy_fy=hpps(A,B,b,i,u,p,fpy_fy);

    % Accumulate the intermediate result
    qq=qq+hpps_up_fpy_fy;

    % Find h'(u)p
    hp_up=hp(A,B,b,i,u,p);

    % Find (f''(y)h'(u)p)^*f(y)
    fpps_yhpup_fy=fpps(i,d,y,hp_up,fy);

    % Find h'(u)^*(f''(y)h'(u)p)^*f(y)
    hpsu_fpps_yhpup_fy=hps(A,B,b,i,u,fpps_yhpup_fy);

    % Accumulate the intermediate result
    qq=qq+hpsu_fpps_yhpup_fy;
end

% Check the difference
if norm(qq-q) > 1e-12
    fprintf('Error when evaluating Newton\n');
end
fprintf('Evaluation of Newton complete.\n');

% Check the Steihaug-Toint process
delta=.3;
max_iter=1000;
eps_cg=1e-8;
s=getStep(A,B,b,d,u,delta,max_iter,eps_cg);

% Set the initial values
s_k=zeros(size(u));
g=getGradient(A,B,b,d,u);
g_k=g;
v_k=g_k;
p_k=-v_k;

for i=1:max_iter
    kappa= GaussNewton(A,B,b,d,u,p_k)*p_k';
    if kappa < 0
	sigma=(-s_k*p_k'+sqrt((s_k*p_k')^2+ ...
	    norm(p_k)^2*(delta^2-norm(s_k)^2)))/ norm(p_k)^2;
	s_k=s_k+sigma*p_k
	break;
    end

    alpha= g_k*v_k'/kappa;

    if norm(s_k+alpha*p_k)>delta
	sigma=(-s_k*p_k'+sqrt((s_k*p_k')^2+ ...
	    norm(p_k)^2*(delta^2-norm(s_k)^2)))/ norm(p_k)^2;
	s_k=s_k+sigma*p_k;
	break;
    end

    s_k = s_k+alpha*p_k;
    old_inner= g_k*v_k';
    g_k = g_k+alpha*GaussNewton(A,B,b,d,u,p_k);

    if norm(g_k)/norm(g) < eps_cg 
	break;
    end

    v_k = g_k;
    beta = g_k*v_k'/old_inner;
    p_k = -v_k + beta * p_k;
end

% Check the difference
if norm(s-s_k) / norm(s_k) > 1e-12
    fprintf('Error when evaluating getStep\n');
end
fprintf('Evaluation of getStep complete.\n');
