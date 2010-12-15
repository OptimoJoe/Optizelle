% Creates the dense BFGS inverse operator given a series of histories of the
% form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function H=build_dense_bfgs_inv(Y,S)

% Get the size of Y 
[m k]=size(Y);

% In this function, we need to reverse the order of Y and S
Y=Y(:,k:-1:1);
S=S(:,k:-1:1);

% Build the "R" matrix
R=zeros(k);
for i=1:k,for j=1:k
    if i<=j
	R(i,j)=S(:,i)'*Y(:,j);
    end
end,end

% Build the "D" matrix
D=zeros(k);
for i=1:k
    D(i,i)=S(:,i)'*Y(:,i);
end

% Build the initial inverse Hessian approximation
H=eye(m);

% Build the whole inverse Hessian approximation
H=H+[S H*Y] * [inv(R)'*(D+Y'*H*Y)*inv(R) -inv(R)';-inv(R) zeros(k)]*[S';Y'*H];
%H=H+[S H*Y] * [R'\((D+Y'*H*Y)*(R\S'))-R'\(Y'*H);-R'\S'];
