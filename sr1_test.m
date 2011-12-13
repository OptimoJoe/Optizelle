% A test of the SR1 and inverse SR1 operators

% Set the random seed
randn('seed',1);

% Set the size of the problem (number of variables)
m=500;

% Set the number of previous points stored  (BFGS uses 1 minus this for the
% history)
n=100;

% Create a list of previous points
X=randn(m,n);
%X=eye(m,n);
%X=[ [1;0;0;0;0] [0;1;0;0;0] ];
%X=[ [1;0;0;0;0] [0;1;0;0;0] [0;0;1;0;0] ];
%p=primes(1000);
%X=reshape(p(1:15),5,3);
%X=reshape(p(1:10),5,2);
	
% Create a list of previous gradients
G=randn(m,n);
%G=eye(m,n);
%G=[ [1;0;0;0;0] [0;1;0;0;0] ];
%G=[ [1;0;0;0;0] [0;1;0;0;0] [0;0;1;0;0] ];
%G=reshape(p(16:30),5,3);
%G=reshape(p(11:20),5,2);

% Determine the difference between the previous points and gradients
S=zeros(m,n-1);
Y=zeros(m,n-1);
for i=n:-1:2
    S(:,n-i+1)=X(:,i)-X(:,i-1);
    Y(:,n-i+1)=G(:,i)-G(:,i-1);
end

% Determine a random direction
p=randn(m,1);

% Apply the SR1 operator
ptild=SR1(Y,S,p);

% Apply the inverse BFGS operator
phat=SR1inv(Y,S,ptild);

% Check positive definateness
alpha=p'*ptild;

% Print some error results
fprintf('The relative error in the operator and inverted direction is: %e\n',norm(p-phat)/norm(p));

% Build matrices that correspond to these operators
H=zeros(m);
Hinv=zeros(m);
for i=1:m
    % Create the ith cannonical vector
    ei=zeros(m,1);
    ei(i)=1;

    % Apply the operators to these cannonical vectors and store the result
    H(:,i)=SR1(Y,S,ei);
    Hinv(:,i)=SR1inv(Y,S,ei);
end

fprintf('In the sparse operators, H*Hinv=%e\n',norm(H*Hinv-eye(m),'fro')/...
    sqrt(m));

% Compute the dense operators directly
H2=build_dense_sr1(Y,S);
Hinv2=build_dense_sr1_inv(Y,S);

% Print out some error information about these operators
fprintf('In the dense operators, H*Hinv=%e\n',norm(H2*Hinv2-eye(m),'fro')/...
    sqrt(m));
