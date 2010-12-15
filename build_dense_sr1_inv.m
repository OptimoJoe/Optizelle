% Creates the dense SR1 inverse operator given a series of histories of the form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function H=build_dense_sr1_inv(Y,S)

% Get the size of Y 
[m k]=size(Y);

% From the first stored history, build the operator
H=eye(m);
for i=k:-1:1
    q=S(:,i)-H*Y(:,i);
    H=H+q*q'/(q'*Y(:,i));
end
