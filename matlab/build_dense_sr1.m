% Creates the dense SR1 operator given a series of histories of the form
% Y=[y_k y_{k-1) ...]
% S=[s_k s_{k-1) ...]

function B=build_dense_sr1(Y,S)

% Get the size of Y 
[m k]=size(Y);

% From the first stored history, build the operator
B=eye(m);
for i=k:-1:1
    q=Y(:,i)-B*S(:,i);
    B=B+q*q'/(q'*S(:,i));
end
