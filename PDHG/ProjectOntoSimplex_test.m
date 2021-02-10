%% Unittests for ProjectOntoSimplex
x = [2,0];
y = ProjectOntoSimplexFast(x,1e-6);
err = y-[1,0]

x = [0,2];
y = ProjectOntoSimplexFast(x,1e-6);
err = y-[0,1]

x = [2,2];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [1/2,1/2]

x = [2,-2];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [1,0]

x = [-2,-2];
y = ProjectOntoSimplexFast(x,1e-6);
err = y-[1/2,1/2]

x = [-2,2];
y = ProjectOntoSimplexFast(x,1e-6);
err = y-[0,1]

x = [0,0];
y = ProjectOntoSimplexFast(x,1e-6);
err = y-[1/2,1/2]

x = [1/4,0];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [5/8,3/8]

x = [0,1/4];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [3/8,5/8];

x = [3/4,0];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [7/8, 1/8]

x = [0,3/4];
y = ProjectOntoSimplexFast(x,1e-6);
err = y - [1/8,7/8]