m = 10000;
n = 1000;

A = randn(m,n); % pretty large matrix
b = randn(m,1); 

% A = [ 1 3 ; 2 -2 ];
% b = [ 1 ; 2 ];

%% true solution

tic
x_solution = A \ b;
toc

%% pdhg

rho2 = normest( A'*A ); % THIS takes more time than solving...
r1 = 0.5;
r2 = 1/rho2/r1;
steps = 70;

toc
x_pdhg = PDHG_solve( A, b, r1, r2, steps );
toc

norm(x_solution - x_pdhg) / norm(x_solution)

sum(x_pdhg)

sum(x_pdhg < 0)

%% conjugate gradients

tic
x_cg = pcg(A'*A, A'*b, 1e-12, steps);
toc

norm(x_solution - x_cg) / norm(x_solution)