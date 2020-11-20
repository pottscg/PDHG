function B_new = update_B(B,A,X,XtX,r1,r2,steps)

    [n,p] = size(B);
    B_new = B;
    [m,n] = size(X);
%     if m >= n
%         rho2 = normest( X'*X ); 
%     else
%         rho2 = normest( X*X' );
%     end
%     r1 = 0.5;
%     r2 = 1/rho2/r1;
%     steps = 10;
%     XtX = X'*X;
    vk = zeros(n,1);
    alpha_k_sqsums = sum(A.*A,2);
    for k = 1:p
        if alpha_k_sqsums(k) > eps
            % update b_k via algorithm, otherwise no update
           B_hat = B_new;
           B_hat(:,k) = [];
           A_hat = A;
           A_hat(k,:) = [];
           vk = (X - X*B_hat*A_hat)*A(k,:)';
           vk = (1/alpha_k_sqsums(k))*sum(vk,2);
           B_new(:,k) = PDHG_solve(X,XtX,vk,r1,r2,steps);
        end
    end
end