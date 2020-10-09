function xn = PDHG_solve(A, AtA, b,r1,r2,steps)
% Solves the min |Ax-b| problem using
% the primal dual hybrid gradient algorithm
% Input: m by n matrix A, m vector b
% Ouput: n vector x


xn = sparse(size(A,2),1);

xbarn = xn;
Atb = A'*b;
% AtA = A'*A;
gradient = AtA*xbarn - Atb;

for k = 1:steps

    xn_last = xn;

    gradient = (gradient + AtA*xbarn - Atb)/(1+r1);
    
    xstar = xn - r1*r2*gradient;
    
    xn = sparse(ProjectOntoSimplex(xstar));
    
    xbarn = 2*xn - xn_last;
    
end

end
