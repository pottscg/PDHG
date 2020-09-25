%% A_NEW = UPDATE_A(A,X,XB)
% Applies PDHG algorithm with simplicial constraints to
% find an update to A for a fixed B in Archetypal Analysis
%
% Inputs: A - reconstruction weights
%         X - data matrix
%         XB - archetypes for a fixed B
% Outputs: A_new - updated reconstruction weights for fixed B

function A_new = update_A(A,X,XB)

    [~,n] = size(X);
    A_new = zeros(size(A));
    rho2 = normest( XB'*XB ); 
    r1 = 0.5;
    r2 = 1/rho2/r1;
    steps = 10;

    for i = 1:n
        A_new(:,i) = PDHG_solve(-XB,-X(:,i),r1,r2,steps);
    end

end








