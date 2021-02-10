profile clear
profile on


N = 1e5;


x1 = 100*randn(N,1)/N;
x1(x1<1) = 0;
x1 = x1 + randn(N,1)/N/4;



tic 
    y1_2 = ProjectOntoSimplexFast(x1,1e-1);
toc

tic 
    y1_1 = ProjectOntoSimplex(x1);
toc

profile off
profile viewer