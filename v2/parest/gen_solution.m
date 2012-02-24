% Returns four functions that correspond to the solution operator,
% its derivative, the adjoint of its derivative, and the adjoint of
% its second derivative after being given the first direction.

function h=gen_solution(A,B,b)

    % Compute the size of the operators Ai and B
    n=size(B,1);

    % Compute the number of Ai operators
    m=length(A);

    % Basic solution operator 
    f=@(i,u)my_f(i,A,B,b,u);

    % Derivative of the solution operator.
    fp=@(i,u,du)my_fp(i,A,B,b,u,du);

    % Piece of the adjoint of the derivative of the solution operator.  
    fps=@(i,u,dy)my_fps(i,A,B,b,u,dy);

    % Adjoint of the second derivative of the solution operator in the
    % direction du applied in the direction dy.  If y does not equal h(u),
    % then this routine is incorrect.
    fpps=@(i,u,du,dy)my_fpps(i,A,B,b,u,du,dy);

    % Build the entire operator
    h.f=f;
    h.fp=fp;
    h.fps=fps;
    h.fpps=fpps;
    h.max_index=length(b);
end

% Computes the operator u1*A1+...+um*Am+B
function C=AuB(A,B,u)
    C=B;
    for i=1:length(u)
    	C=C+A{i}*u(i);
    end
end

% Create the operation A(.)y 
function z=Adoty(A,y,delta_u)
    z=zeros(size(y));	
    for i=1:length(u)
    	z=z+u(i)*A*y{i}
    end
end

% Define the operator [A(.)y]*
function v=Adoty_adj(A,y,delta_y)
    v=zeros(length(A),1);
    for i=1:length(v)
	v(i)=(A{i}*delta_y)'*y;
    end
end

% This solves the governing equations.  It's very important that this function
% caches results effectively
function y=my_f(idx,A,B,b,u)
    % These are two cell arrays that contain saved states and the points that
    % generated them
    persistent saved_y
    persistent saved_u

    % Initialize the first time we run
    if isempty(saved_y)
    	saved_y=cell(0);
	saved_u=cell(0);
    end

    % Check if we have already calculated the state 
    mynorm=@(x)sqrt(innr(x,x));
    if idx <= length(saved_y)
    	if mynorm(saved_u{idx}-u) < 1e-16
	    y=saved_y{idx};
	    return
	end
    end

    % Generate the operator
    C=AuB(A,B,u);

    % Find the state
    y=C\b{idx};

    % Save the state
    saved_y{idx}=y;
    saved_u{idx}=u;
end

% Finds the derivative of the solution operator
function dy=my_fp(idx,A,B,b,u,du)

    % Find the state solution
    y=my_f(idx,A,B,b,u);

    % Find A(du)y
    z=zeros(size(y));
    for i=1:length(u)
	z=z+du(i)*(A{i}*y);
    end
   
    % Generate the full governing equations
    C=AuB(A,B,u);

    % Find dy
    dy=-C\z;
end

% Finds the adjoint of the derivative of the solution operator
function du=my_fps(idx,A,B,b,u,dy)

    % Find the state solution
    y=my_f(idx,A,B,b,u);

    % Find the adjoint solution 
    C=AuB(A,B,u)';
    z=C\dy;

    % Find -[A(.)y]*z
    du=zeros(size(u));
    for i=1:length(du)
    	du(i)=-(A{i}*y)'*z;
    end
end

% Finds the second derivative applied to du, adjoint, applied to dy
function dv=my_fpps(idx,A,B,b,u,du,dy)

    % Find the adjoint solution
    C=AuB(A,B,u)';
    z=C\dy;

    % Find A(du)*z
    A_du_z=zeros(size(z));
    for i=1:length(u)
	A_du_z=A_du_z+du(i)*(A{i}'*z);
    end

    % Find h'(u)*A_du_z
    hps_A_du_z=my_fps(idx,A,B,b,u,A_du_z);

    % Find h'(u)du
    hp_du=my_fp(idx,A,B,b,u,du);

    % Find [A(.)h'(u)du]*z
    A_hp_du_z=zeros(size(u));
    for i=1:length(u)
    	A_hp_du_z(i)=(A{i}*hp_du)'*z;
    end

    % Combine everything above for the final adjoint solution
    dv=-A_hp_du_z-hps_A_du_z;
end
