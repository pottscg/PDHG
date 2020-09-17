function xn = PDHG_solve(A,b,r1,r2,steps)
% Solves the min |Ax-b| problem using
% the primal dual hybrid gradient algorithm
% Input: m by n matrix A, m vector b
% Ouput: n vector x


xn = zeros(size(A,2),1);

xbarn = xn;
Atb = A'*b;
AtA = A'*A;
gradient = AtA*xbarn - Atb;

for k = 1:steps

    xn_last = xn;

    gradient = (gradient + AtA*xbarn - Atb)/(1+r1);
    
    xn = xn - r1*r2*gradient;
    
    xn = ProjectOntoSimplex(xn);
    
    xbarn = 2*xn - xn_last;
    
end

end
