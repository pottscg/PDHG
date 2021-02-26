addpath('PCHA/');

filename = 'cbcl_FJLT_10MorupAA_unnormalizedGaussian';

load('DATASETS/cbcl_FJLT.mat');
% load(strcat('RESULTS/MORUP/',filename,'.mat'));
[m,n] = size(X);
noc=10; % Number of archetypes

delta=0;
opts.maxiter=1000;
opts.conv_crit=1e-6;
U = 1:n;
I = 1:n;

[XC,S,C,SSE,varexpl]=PCHA(X,noc,I,U,delta,opts);

L2_Model_err = zeros(85,20);
L2_MORUPRECON_ERR = zeros(85,20);
L2_ORIGINAL_ERR = zeros(85,20);
L2_MORUPC_ERR = zeros(85,20);
L2_MORUPS_ERR = zeros(85,20);
L2_MORUPXC_ERR = zeros(85,20);
Recon_Morup = X*C*S;

for seed = 1:20
    rng(seed);
    for k = 1:85
        d = 3*k;
        % d = 25; % Number of dimensions to project onto

        %Initialize sketches


%         %create d x m projection matrix 
% %         PSI = reshape(rand(d*m,1),[d,m]);
%         PSI = randn(d,m)
%         
%         %normalize over rows
% %         for j = 1:d
% %             PSI(j,:) = PSI(j,:)./sum(PSI(j,:));
% %         end

        PSI = Sketch_Projection(m,d,'Gaussian');

        %project data
        Y = zeros(d,n);

        for i = 1:n
            y = PSI*X(:,i);
            Y(:,i) = reshape(y,[d,1]);
        end



        U=1:size(Y,2); % Entries in Y used that is modelled by the AA model
        I=1:size(Y,2); % Entries in Y used to define archetypes
        % if two expensive to useall entries for I find N relevant observations by
        % the following procedure:
        % N=100;
        % I=FurthestSum(Y,N,ceil(rand*size(Y,2)));

        %Initialize sketches

%         disp(strcat('-------- seed val - ',num2str(seed),'--- dim - ',num2str(d)));
        
        % Use PCHA.m
        [YC,Stilda,Ctilda,SSE,varexpl]=PCHA(Y,noc,I,U,delta,opts);

        % Calculate archetypes from full data
        XC_sketches = X*Ctilda;
        Recon_sketches = XC_sketches*Stilda;
        
%         % Rearrange based on archetype match
%         XC_sketches_arranged = XC_sketches;
%         Ctilda_arranged = Ctilda;
%         Stilda_arranged = Stilda;
%         for j = 1:noc
%             Archetype = repmat(XC_sketches(:,j),1,noc);
%             [M,I] = min(vecnorm(Archetype - XC,2));
%             XC_sketches_arranged(:,I) = XC_sketches(:,j);
%             Ctilda_arranged(:,I) = Ctilda(:,j);
%             Stilda_arranged(I,:) = Stilda(j,:);
%         end
%         Recon_sketches = XC_sketches_arranged*Stilda_arranged;

        % Generate L2 errors
        L2_MORUPRECON_ERR(k,seed) = norm(Recon_sketches - Recon_Morup,'fro')/size(Recon_sketches,2);
        L2_ORIGINAL_ERR(k,seed) = norm(Recon_sketches - X,'fro')/size(Recon_sketches,2);
        L2_MORUPC_ERR(k,seed) = norm(Ctilda - C,'fro')/size(Ctilda,2);
        L2_MORUPS_ERR(k,seed) = norm(Stilda - S,'fro')/size(Stilda,2);
        L2_MORUPXC_ERR(k,seed) = norm(XC_sketches- XC, 'fro')/size(XC_sketches,2);
        L2_Model_err(k,seed) = norm(X - Recon_sketches, 'fro')^2;
    end
end
save(strcat('RESULTS/SKETCH_TYPES/',filename,'.mat'),'X','L2_MORUPRECON_ERR','L2_ORIGINAL_ERR',...
        'L2_MORUPC_ERR', 'L2_MORUPS_ERR','L2_MORUPXC_ERR', 'L2_Model_err', 'Recon_Morup');