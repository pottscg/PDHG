function w = ProjectOntoSimplex(v)
% PROJECTONTOSIMPLEX Projets point onto simplex of specifieid radius. 
%
% w = ProjectOntoSimplex(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%   
%   min  ||w-v||_2
%   s.t. sum(w) <=b, w >=0
%
%  That is, performs Euclidean projection of v to the positive simplex of
%   radius b. 
%
% Auher: John Duchi (jduchi@cs.berkley.edu)

v = (v>0).*v;
u = sort(v,'descend');
sv = cumsum(u);
rho = find(u > (sv-1) ./ (1:length(u))',1,'last');
theta = max(0, (sv(rho)-1)/rho);
w = max(v-theta,0);