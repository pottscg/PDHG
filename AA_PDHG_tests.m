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
[S,SSE,muS,SSt]=Supdate(S,XCtX,CtXtXC,muS,SST,SSE,10);
toc;