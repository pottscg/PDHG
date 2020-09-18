function w = ProjectOntoSimplex(v)
% PROJECTONTOSIMPLEX Projects point v onto the probability simplex.

v = v - min(v(:));   % make sure we are in the positive orthant

sum_v = sum(v(:));  % compute the current l1-norm

if sum_v < 1  
  w = v + (1-sum_v)/numel(v);   % if deficient, just uniformly add
else  
  u = sort(v,'descend');      % if excessive, "just shrink"
  sv = cumsum(u);
  rho = find(u > (sv-1) ./ (1:numel(u))',1,'last');
  theta = max(0, (sv(rho)-1)/rho);
  w = max(v-theta,0);
end
