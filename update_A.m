%% A_NEW = UPDATE_A(A,X,XB)
% Applies PDHG algorithm with simplicial constraints to
% find an update to A for a fixed B in Archetypal Analysis
%
% Inputs: A - reconstruction weights
%         X - data matrix
%         XB - archetypes for a fixed B
% Outputs: A_new - updated reconstruction weights for fixed B

function [A_new,SSE] = update_A(A,X,XB,SST, noc)

    [~,n] = size(X);
    e=ones(noc,1);
    A_new = zeros(size(A));
    XBtXB = XB'*XB;
    rho2 = normest( XBtXB ); 
    r1 = 0.0005;
    r2 = 1/rho2/r1;
    steps = 5;
    for i = 1:n
        A_new(:,i) = PDHG_solve(-XB,XBtXB,-X(:,i),A(:,i),r1,r2,steps);
    end
    % morup's convexity update
%     A_new(A_new<0)=0;
%     A_new=A_new./(e*sum(A_new,1));
    SSE=SST-2*sum(sum((XB'*X).*A_new))+sum(sum(XBtXB.*(A_new*A_new'))); 
end








