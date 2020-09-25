profile on

addpath '/home/potts/Archetypal_Analysis/PDHG'

load('/home/potts/Archetypal_Analysis/RESULTS/MORUP/mnist_20MorupAA.mat');

%% Build reconstruction data off Morup's results to compare with additional update

%Wiggle S away from optimal
r = -0.25 + 0.5*rand(size(S));
S = S + r;
X_recon = XC*S;
Original_err = norm(X-X_recon)/norm(X)

%% Test update of A given B fixed

tic;

A_PDHGupdate = update_Agpu(S,X,XC);

toc;

X_PDHG_recon = XC*A_PDHGupdate;
PDHG_err = norm(X-X_PDHG_recon)/norm(X)

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

X_Morup_recon = XC*Supdate;
Morup_err = norm(X-X_Morup_recon)/norm(X)

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

    parfor i = 1:n
        A_new(:,i) = lsqlin(-XB,-X(:,i)',ones(1,size(A,1)), 1, [], [], zeros(size(A,1),1), []);
    end
 
  toc;
 
X_lsqlin_recon = XC*A_new;
Lsqlin_err = norm(X-X_lsqlin_recon)/norm(X)

profile off
profile viewer
