% Creates the dense BFGS operator given a series of histories of the form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function B=build_dense_bfgs(Y,S)

% Get the size of Y 
[m k]=size(Y);

% From the first stored history, build the operator
B=eye(m);
for i=k:-1:1
    B=B-(B*S(:,i)*S(:,i)'*B)/(S(:,i)'*B*S(:,i));
    B=B+Y(:,i)*Y(:,i)'/(Y(:,i)'*S(:,i));
end
