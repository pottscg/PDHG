load('DATASETS/cbcl_FJLT.mat');
filename = 'cbcl_FJLT_Bad_epsp25';
sketch_type = 'Bad';
epsilon = 0.25;

[m,n] = size(X);

JL_Satisfied = zeros(85,20);

for seed = 1:20
    rng(seed);
    for k = 1:85
        d = 3*k;
        

        PSI = Sketch_Projection(m,d,sketch_type);

        %project data
        Y = zeros(d,n);

        for i = 1:n
            y = PSI*X(:,i);
            Y(:,i) = reshape(y,[d,1]);
        end

        result = JL_satisfied(X,Y,epsilon);
        
        JL_Satisfied(k,seed) = any(result == 0, 'all');
        disp(strcat(num2str(d),'_',num2str(seed), '_' ,num2str(JL_Satisfied(k,seed))));
    end
end
save(strcat('RESULTS/JL_Satisfied/',filename,'.mat'),'JL_Satisfied');