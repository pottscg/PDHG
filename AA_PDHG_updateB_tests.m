clear all; 

% profile on

addpath '/home/potts/Archetypal_Analysis/PDHG'

load('/home/potts/Archetypal_Analysis/RESULTS/MORUP/mnist_20MorupAA.mat');

%% Build reconstruction data off Morup's results to compare with additional update

%Wiggle C away from optimal
r = -0.25 + 0.5*rand(size(C));
C = C + r;
XC = X*C;
X_recon = XC*S;
Original_err = norm(X-X_recon)/norm(X)

%% Test update of B given A fixed

tic;

B_PDHGupdate = update_B(C,S,X);

toc;

X_PDHG_recon = X*B_PDHGupdate*S;
PDHG_err = norm(X-X_PDHG_recon)/norm(X)

%% Test C (B) update given S (A) fixed via Morup's method

tic;
%pre-lim set-up
U=1:size(X,2);
I=1:size(X,2);
XSt=X(:,U)*S';
SSt=S*S';
delta=0;
muC=1;
mualpha=1;
SST=sum(sum(X(:,U).*X(:,U)));

XCtX=XC'*X(:,U);
CtXtXC=XC'*XC;

SSE=SST-2*sum(sum(XCtX.*S))+sum(sum(CtXtXC.*SSt));
[Cupdate,SSE,muC,mualpha,CtXtXC,XCupdate]=Cupdate(X(:,I),XSt,XC,SSt,C,delta,muC,mualpha,SST,SSE,10);
toc;

X_Morup_recon = X*Cupdate*S;
Morup_err = norm(X-X_Morup_recon)/norm(X)


% profile off
% profile viewer % shows the html profile
%or
% profsave(profile('info'),'update_Agpu_test_profileresults')% saves an html in current working directory
