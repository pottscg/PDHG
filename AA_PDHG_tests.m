addpath '/home/potts/Archetypal_Analysis/PDHG'

load('/home/potts/Archetypal_Analysis/RESULTS/MORUP/mnist_20MorupAA.mat');

%% Test update of A given B fixed

tic;

A_PDHGupdate = update_A(S,X,XC);

toc;

%% Test S (A) update given C (B) fixed via Morup's method

tic;
%pre-lim set-up
SST=sum(sum(X(:,U).*X(:,U)));
CtXtXC = XC'*XC;
XCtX=XC'*X(:,U);
SSt=S*S'; 
XSt = X*S';
SSE=SST-2*sum(sum(XC.*XSt))+sum(sum(CtXtXC.*SSt));
muS = 1;
%update
[Supdate,SSEupdate,muSupdate,SStupdate]=Supdate(S,XCtX,CtXtXC,muS,SST,SSE,10);
toc;

%% Test A update using Matlab's lsqlin

tic;
A = S;
XB = XC;
[m,n] = size(X);
    A_new = zeros(size(A));
    tic;
    rho2 = normest( XB'*XB ); 
    toc;
    r1 = 0.5;
    r2 = 1/rho2/r1;
    steps = 10;

    for i = 1:n
        A_new(:,i) = lsqlin(-XB,-X(:,i)',ones(1,size(A,1)), 1, [], [], zeros(size(A,1),1), []);
    end
 
  toc;