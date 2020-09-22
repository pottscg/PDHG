%% Unittests for ProjectOntoSimplex
x = [2,0]
y = ProjectOntoSimplex(x)
err = y-[1,0]

x = [0,2]
y = ProjectOntoSimplex(x)
err = y-[0,1]

x = [2,2]
y = ProjectOntoSimplex(x)
err = y - [1/2,1/2]

x = [2,-2]
y = ProjectOntoSimplex(x)
err = y - [1,0]

x = [-2,-2]
y = ProjectOntoSimplex(x)
err = y-[1/2,1/2]

x = [-2,2]
y = ProjectOntoSimplex(x)
err = y-[0,1]

x = [0,0]
y = ProjectOntoSimplex(x)
err = y-[1/2,1/2]

x = [1/4,0]
y = ProjectOntoSimplex(x)

x = [0,1/4]
y = ProjectOntoSimplex(x)

x = [3/4,0]
y = ProjectOntoSimplex(x)

x = [0,3/4]
y = ProjectOntoSimplex(x)