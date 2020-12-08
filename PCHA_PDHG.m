function [XB,A,B,SSE,varexpl]=PCHA_PDHG(X,noc,I,U,delta,varargin)
% Principal Convex Hull Analysis (PCHA) / Archetypal Analysis
%
% Written by Morten M�rup
%
% Usage:
%   [XC,S,C,SSE,varexpl]=PCHA(X,noc,W,I,U,delta,varargin)
%
%   Solves the following PCH/AA problem
%   \|X(:,U)-X(:,I)CS\|_F^2 s.t. |s_j|_1=1, 1-delta<=|c_j|_1<=1+delta,
%   S>=0 and C>=0
%
%
% Input:
% X             data array 
% noc           number of components
% I             Entries of X to use for dictionary in C (default: I=1:size(X,2))
% U             Entries of X to model in S              (default: U=1:size(X,2))
% delta         relaxation of C, i.e. 1-delta<=|C_j|_1<=1+delta (default: delta=0, i.e. |C_j|_1=1)
%
% opts.         Struct containing:
%       C            initial solution (optional) (see also output)
%       S            initial solution (optional) (see also output)
%       maxiter      maximum number of iterations (default: 500 iterations)
%       conv_crit    The convergence criteria (default: 10^-6 relative change in SSE)
%
% Output:
% XC            I x noc feature matrix (i.e. XC=X(:,I)*C forming the archetypes) 
% S             noc x length(U) matrix, S>=0 |S_j|_1=1
% C             length(I) x noc matrix, C>=0 1-delta<=|C_j|_1<=1+delta
% SSE           Sum of Squares Error
% varexpl       Percent variation explained by the model
%
% Copyright (C) Morten M�rup and Technical University of Denmark, 2010
%
% Updated 8/10-2018 to reflect in the description that this code does not
% support missing data (for missing data see PCHAsparse.m)
% Updated 9/10-2018 to handle when noc is specified as 1 (sum(*) changed to sum(*,1) in code)

warning('off','MATLAB:dispatcher:InexactMatch')
if nargin>=6, opts = varargin{1}; else opts = struct; end
conv_crit=mgetopt(opts,'conv_crit',10^-6);
maxiter=mgetopt(opts,'maxiter',500);

if nargin<5
    delta=0;
end
if nargin<4
    U=1:size(X,2);
end
if nargin<3
    I=1:size(X,2);
end

%sum square total
SST=sum(sum(X(:,U).*X(:,U)));

% Initilize B 
if isfield(opts,'B')
    B = opts.B;    
else
   % Initialize by furthest sum   
   i=FurthestSum(X(:,I),noc,ceil(length(I)*rand));
   B=sparse(i,1:noc,ones(1,noc),length(I),noc);    
end
% archetypes
XB=X(:,I)*B; 

% steps
muA=1;
muB=1;
mualpha=1;

% Initilize A 
if isfield(opts,'A')
    A=opts.A;    
    BtXtXB=XB'*XB;    
    XAt=X(:,U)*A';
    AAt=A*A';           
    SSE=SST-2*sum(sum(XB.*XAt))+sum(sum(BtXtXB.*AAt));                
else   
%     XBtX=XB'*X(:,U);
%     BtXtXB=XB'*XB;    
    A=-log(rand(noc,length(U)));
    A=A./(ones(noc,1)*sum(A,1));    
%     AAt=A*A';
%     SSE=SST-2*sum(sum(XBtX.*A))+sum(sum(BtXtXB.*AAt));                
%     [A,SSE,muA,AAt]=Supdate(A,XBtX,BtXtXB,muA,SST,SSE,25);
    [A,SSE] = update_A(A,X(:,U),XB,SST);
end

% looking at initial archetypes and weights
figure; 
hold on;
scatter(X(1,:),X(2,:),2,'filled'); 
xlim([0 3]); ylim([0 3]);
scatter(XB(1,:),XB(2,:),25,[0.85, 0.325, 0.098],'filled');
hold off;


% Set PCHA parameters
iter=0;
dSSE=inf;
t1=cputime;
varexpl=(SST-SSE)/SST;

% B update parameters - fixed for entire run
    XtX = X'*X;
    [m,n] = size(X);
    if m >= n
        rho2 = normest( XtX ); 
    else
        rho2 = normest( X*X' );
    end
    r1 = 0.5;
    r2 = 1/rho2/r1;
    steps = 10;

% Display algorithm profile
disp([' '])
disp(['Principal Convex Hull Analysis / Archetypal Analysis'])
disp(['A ' num2str(noc) ' component model will be fitted']);
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s','Iteration','Expl. var.','Cost func.','Delta SSEf.','muC','mualpha','muS',' Time(s)   ');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+');

told = 0; %CP - added to avoid stoprun error
t1 = 0;
while abs(dSSE)>=conv_crit*abs(SSE) && iter<maxiter && varexpl<0.9999
    if mod(iter,100)==0
         disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    SSE_old=SSE;
        
    % B (and alpha) update
%     XAt=X(:,U)*A';
%     [B,SSE,muB,mualpha,BtXtXB,XB]=Cupdate(X(:,I),XAt,XB,AAt,B,delta,muB,mualpha,SST,SSE,10);

    B = update_B(B,A,X,XtX,r1,r2,steps);
    XB = X*B;
    % A update    
%     XBtX=XB'*X(:,U);    
%     [A,SSE,muA,AAt]=Supdate(A,XBtX,BtXtXB,muA,SST,SSE,10);  
    [A,SSE] = update_A(A,X,XB,SST);
  
%     SSE=SST-2*sum(sum((XB'*X).*A))+sum(sum(BtXtXB.*AAt)); 
    % Evaluate and display iteration
    dSSE=SSE_old-SSE;
    t1=cputime;
    if rem(iter,1)==0  
        pause(0.000001);
        varexpl=(SST-SSE)/SST;
        fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,full(varexpl),full(SSE),full(dSSE/abs(SSE)),muB,mualpha,muA,t1);
    end
end

% display final iteration
varexpl=(SST-SSE)/SST;
disp(dline);
disp(dline);
fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,full(varexpl),full(SSE),full(dSSE/abs(SSE)),muB,mualpha,muA,t1);

% sort components according to importance
[val,ind]=sort(sum(A,2),'descend');
A=A(ind,:);
B=B(:,ind);
XB=XB(:,ind);

% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end


% % -------------------------------------------------------------------------
% function [S,SSE,muS,SSt]=Supdate(S,XCtX,CtXtXC,muS,SST,SSE,niter)
%     
%     [noc,J]=size(S);
%     e=ones(noc,1);
%     for k=1:niter
%         SSE_old=SSE;
%         g=(CtXtXC*S-XCtX)/(SST/J);         
%         g=g-e*sum(g.*S,1);
%         stop=0;
%         Sold=S;
%         while ~stop
%             S=Sold-g*muS;
%             S(S<0)=0;
%             S=S./(e*sum(S,1));            
%             SSt=S*S';          
%             SSE=SST-2*sum(sum(XCtX.*S))+sum(sum(CtXtXC.*SSt));            
%             if SSE<=SSE_old*(1+1e-9)
%                 muS=muS*1.2;
%                 stop=1;
%             else
%                 muS=muS/2;
%             end
%         end
%     end
% 
% %--------------------------------------------------------------------
% function [C,SSE,muC,mualpha,CtXtXC,XC]=Cupdate(X,XSt,XC,SSt,C,delta,muC,mualpha,SST,SSE,niter)
%                                        
%     [J,noc]=size(C);
%     if nargin<12
%         niter=1;
%     end   
%     if delta~=0
%         alphaC=sum(C);    
%         C=C*diag(1./alphaC);
%     end
%     e=ones(J,1);
%     XtXSt=X'*XSt;
%     for k=1:niter
%         
%         % Update C        
%         SSE_old=SSE;        
%         g=(X'*(XC*SSt)-XtXSt)/SST;       
%         if delta~=0
%             g=g*diag(alphaC);
%         end
%         g=g-e*sum(g.*C,1);        
%         stop=0;
%         Cold=C;
%         while ~stop
%             C=Cold-muC*g;
%             C(C<0)=0;
%             
%             nC=sum(C,1)+eps;            
%             C=C*sparse(1:noc,1:noc,1./nC);
%             if delta~=0
%                 Ct=sparse(C*diag(alphaC));
%             else
%                 Ct=sparse(C);
%             end
%             XC=X*Ct;
%             CtXtXC=XC'*XC;
%             SSE=SST-2*sum(sum(XC.*XSt))+sum(sum(CtXtXC.*SSt));                                    
%             if SSE<=SSE_old*(1+1e-9)
%                 muC=muC*1.2;
%                 stop=1;
%             else
%                 muC=muC/2;
%             end
%         end        
%         
%         % Update alphaC        
%         SSE_old=SSE;
%         if delta~=0                                                           
%             g=(diag(CtXtXC*SSt)'./alphaC-sum(C.*XtXSt,1))/(SST*J);                       
%             stop=0;
%             alphaCold=alphaC;
%             while ~stop
%                 alphaC=alphaCold-mualpha*g;
%                 alphaC(alphaC<1-delta)=1-delta;
%                 alphaC(alphaC>1+delta)=1+delta;                            
%                 XCt=XC*diag(alphaC./alphaCold);
%                 CtXtXC=XCt'*XCt;            
%                 SSE=SST-2*sum(sum(XCt.*XSt))+sum(sum(CtXtXC.*SSt));                
%                 if SSE<=SSE_old*(1+1e-9)  
%                     mualpha=mualpha*1.2;
%                     stop=1;
%                     XC=XCt;
%                 else
%                     mualpha=mualpha/2;
%                 end
%             end  
%         end
%     end
%     if delta~=0
%         C=C*diag(alphaC);
%     end
% 
