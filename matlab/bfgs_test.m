% A test of the BFGS and inverse BFGS operators

% Set the random seed
%randn('seed',1);

% Set the size of the problem (number of variables)
m=50;

% Set the number of previous points stored  (BFGS uses 1 minus this for the
% history)
n=10;

% Create a list of previous points
%X=randn(m,n);
%X=eye(m,n);
%X=[ [1;0;0;0;0] [0;1;0;0;0] ];
%X=[ [1;0;0;0;0] [0;1;0;0;0] [0;0;1;0;0] ];
%p=primes(1000);
%X=reshape(p(1:15),5,3);
%X=reshape(p(1:10),5,2);
	
% Create a list of previous gradients
%G=randn(m,n);
%G=eye(m,n);
%G=[ [1;0;0;0;0] [0;1;0;0;0] ];
%G=[ [1;0;0;0;0] [0;1;0;0;0] [0;0;1;0;0] ];
%G=reshape(p(16:30),5,3);
%G=reshape(p(11:20),5,2);

% Be more careful about generating a list of previous points
X=[randn(m,1) zeros(m,n-1)];
G=[randn(m,1) zeros(m,n-1)];
for j=2:n
    X(:,j)=randn(m,1);
    s=X(:,j)-X(:,j-1);
    G(:,j)=randn(m,1);
    y=G(:,j)-G(:,j-1);
    while(s'*y<=0)
	G(:,j)=randn(m,1);
	y=G(:,j)-G(:,j-1);
    end
end

% Determine the difference between the previous points and gradients
S=zeros(m,n-1);
Y=zeros(m,n-1);
for i=n:-1:2
    S(:,n-i+1)=X(:,i)-X(:,i-1);
    Y(:,n-i+1)=G(:,i)-G(:,i-1);
end

% Determine a random direction
p=randn(m,1);

% Apply the BFGS operator
ptild=BFGS(Y,S,p);

% Apply the inverse BFGS operator
phat=BFGSinv(Y,S,ptild);

% Check positive definateness
alpha=p'*ptild;

% Print some error results
fprintf('The relative error in the operator and inverted direction is: %e\n',norm(p-phat)/norm(p));
fprintf('The H-inner product of p is: %e\n',alpha);

% Build matrices that correspond to these operators
H=zeros(m);
Hinv=zeros(m);
for i=1:m
    % Create the ith cannonical vector
    ei=zeros(m,1);
    ei(i)=1;

    % Apply the operators to these cannonical vectors and store the result
    H(:,i)=BFGS(Y,S,ei);
    Hinv(:,i)=BFGSinv(Y,S,ei);
end

fprintf('In the sparse operators, H*Hinv=%e\n',norm(H*Hinv-eye(m),'fro')/...
    sqrt(m));

% Assume that we have a history of two vectors.  Build the inverse Hessian
% operator directly
%s=S(:,1);
%y=Y(:,1);
%rho=1/(s'*y);
%Hinv2=(eye(m)-rho*s*y')*(eye(m)-rho*y*s')+rho*s*s';

% Do the same for the forward operator
%H2=eye(m)-s*s'/(s'*s) + y*y'/(y'*s);

% Compute the dense operators directly
H2=build_dense_bfgs(Y,S);
Hinv2=build_dense_bfgs_inv(Y,S);

% Print out some error information about these operators
fprintf('In the dense operators, H*Hinv=%e\n',norm(H2*Hinv2-eye(m),'fro')/...
    sqrt(m));

% Compute the dense operators using the compact representation
H3=build_compact_bfgs(Y,S);
Hinv3=build_compact_bfgs_inv(Y,S);

% Print out some error information about these operators
fprintf('In the compact operators, H*Hinv=%e\n',norm(H3*Hinv3-eye(m),'fro')/...
    sqrt(m));
