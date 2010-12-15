% A test of the BFGS and inverse BFGS operators

% Set the random seed
%randn('seed',10);

% Set the size of the problem (number of variables)
m=200;

% Set the number of previous points stored  (BFGS uses 1 minus this for the
% history)
n=100;

% Be more careful about generating a list of previous points
%X=[randn(m,1) zeros(m,n-1)];
X=zeros(m,n);X(1)=1;
%G=[randn(m,1) zeros(m,n-1)];
G=zeros(m,n);G(1)=1;
eval_history=zeros(m,n-1);
for j=2:n
    %X(:,j)=randn(m,1);
    X(j,j)=randn;
    s=X(:,j)-X(:,j-1);
    %G(:,j)=randn(m,1);
    G(j,j)=randn;
    y=G(:,j)-G(:,j-1);
    while(s'*y<=0)
	%G(:,j)=randn(m,1);
	G(j,j)=randn;
	y=G(:,j)-G(:,j-1);
    end

    % Determine the difference between the previous points and gradients
    S=zeros(m,j-1);
    Y=zeros(m,j-1);
    for i=j:-1:2
	S(:,j-i+1)=X(:,i)-X(:,i-1);
	Y(:,j-i+1)=G(:,i)-G(:,i-1);
    end

    % Once we have these histories, find the eigenvalues of the forward operator
    H=build_dense_bfgs(Y,S);
    [v d]=eig(H);
    eval_history(:,j-1)=real(diag(d));
end

clf
for i=1:n-1
    semilogy(ones(m,1)*i,eval_history(:,i),'x')
    hold on
end
hold off
xlabel('Number of stored vectors')
ylabel('Eigenvalues')
title(sprintf('Eigenvalues of BFGS operator of size: %d',m));

Hinv=build_dense_bfgs_inv(Y,S);
norm(H*Hinv-eye(m))
