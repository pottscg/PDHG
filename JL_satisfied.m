% Calculates the distances between pairwise points for sketched data 

function result = JL_satisfied(X,S, epsilon)

    [m,n] = size(X);
%     
%     pairwise_X_dist = zeros(n,n);
%     pairwise_S_dist = zeros(n,n);
    result = zeros(n,n);
    
    for j = 1:n
        for k = j:n
            X_dist = norm(X(:,j)-X(:,k))^2;
            S_dist = norm(S(:,j)-S(:,k))^2;
            if (1+epsilon)*X_dist < S_dist && (1-epsilon)*X_dist > S_dist
                result(j,k) = 1;
                result(k,j) = 1;
            end
            if j==k
                result(j,k) = 1;
            end
        end
    end

    
end