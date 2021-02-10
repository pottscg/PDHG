function w = ProjectOntoSimplexFast(v,tol)
% PROJECTONTOSIMPLEXFAST Projects point v onto the probability simplex.

N = numel(v);
min_v = min(v(:));
if min_v < 0
    v = v - min_v;   % make sure we are in the positive orthant
end
sum_v = sum(v(:));  % compute the current l1-norm

if sum_v < 1
    w = v + (1-sum_v)/N;   % if deficient, just uniformly add
else
    theta_min = (sum_v - 1) / N; % same amount shrunk from all components
    theta_max = 100*theta_min; %(sum_v - 1); % amount shrunk from a single nnz component
    theta = (theta_max + theta_min)/2;
    %     w = max(v-theta,zeros(size(v)));
    w = v - min(v,theta);
    sum_v = sum(w,"all");
    while abs(sum_v - 1) > tol
%         disp(theta_min + " < " + theta + " < " + theta_max);
        if sum_v > 1
            theta_min = theta;
            theta = (theta + theta_max)/2;
        else
            theta_max = theta;
            theta = (theta + theta_min)/2;
        end
        w = v - min(v,theta);
        sum_v = sum(w,"all");
        
    end
    %     u = sort(v(:),'descend');      % if excessive, "just shrink"
    %     sv = cumsum(u);
    %     rho = find(u > (sv-1) ./ (1:numel(u))',1,'last');
    %     theta = max(0, (sv(rho)-1)/rho);
    %     w = max(v-theta,zeros(size(v)));
end
