% Creates the dense BFGS operator given a series of histories of the
% form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function B=build_dense_bfgs_inv(Y,S)

% Get the size of Y 
[m k]=size(Y);

% In this function, we need to reverse the order of Y and S
Y=Y(:,k:-1:1);
S=S(:,k:-1:1);

% Build the "L" matrix
L=zeros(k);
for i=1:k,for j=1:k
    if i>j
	L(i,j)=S(:,i)'*Y(:,j);
    end
end,end

% Build the "D" matrix
D=zeros(k);
for i=1:k
    D(i,i)=S(:,i)'*Y(:,i);
end

% Build the initial Hessian approximation
B=eye(m);

% Build the whole Hessian approximation
B=B-[B*S Y]*([S'*B*S L;L' -D]\[S'*B;Y']);
