% Solves the min |Ax-b| problem using
% the primal dual hybrid gradient algorithm
% Input: m by n matrix A, m vector b
% Ouput: m vector x

function xn = PDHG_solve(A,b,r1,r2,theta)

steps = 15;
xn = zeros(size(A,2),1);
yn = zeros(size(A,1),1);
xbarn = xn;
Atb = A'*b;
AtA = A'*A;
gradient = zeros(size(A,2),steps+1);
gradient(:,1) = AtA*xbarn - Atb;

for k = 1:steps
%     yn = (1/(r1+1))*yn + (r1/(r1+1))*(A*xbarn - b);
%     xn_last = xn;
%     xn = xn - r2*A'*yn;
%     xbarn = xn +theta*(xn - xn_last);
%     
% without dual variable
    
    xn_last = xn;
    % calculate gradient step for current k time step
    gradient_step = zeros(size(A,2),1);
    for j = 0:k
        gradient_step =  gradient_step + gradient(:,j+1)*(1/(1+r1)^(k-j+1));
    end
    xn = xn - r1*r2*gradient_step;
    
    %projection todo
    
    xbarn = (1+theta)*xn - theta*xn_last;
    gradient(:,k+1) = AtA*xbarn - Atb;
    
end

end