% Creates the dense BFGS inverse operator given a series of histories of the
% form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function H=build_dense_bfgs_inv(Y,S)

% Get the size of Y 
[m k]=size(Y);

% From the first stored history, build the operator
H=eye(m);
for i=k:-1:1
    rho=1/(Y(:,i)'*S(:,i));
    H=(eye(m)-rho*S(:,i)*Y(:,i)')*H*(eye(m)-rho*Y(:,i)*S(:,i)');
    H=H+rho*S(:,i)*S(:,i)';
end
