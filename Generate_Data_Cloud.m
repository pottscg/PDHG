%% 3 clusters
dim_data = 100;
num_datapts = 900;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate multidimensional mean values
mus = mvnrnd(zeros(dim_data,1)',eye(dim_data), 3);
muA = mus(1,:);
muA(1:15) = 0;
muA = muA./mean(muA);
muA(abs(muA) > 50) = 0;

muB = mus(2,:);
muB(16:30) = 0;
muB = muB./mean(muB);
muB(abs(muB) > 40) = 0;

muC = mus(3,:);
muC(50:65) = 0;
muC = muC./mean(muC);
muC(abs(muC) > 30 ) = 0;

% generate multidimensional covariance matrix
sigma = diag(mvnrnd(zeros(dim_data,1),eye(dim_data),1));
sigma = abs(sigma);  %enforce symmetric positive definite
% sigma = sigma./mean(sigma(:));

X(:,1:300) = mvnrnd(muA,sigma,300)';
X(:,301:600) = mvnrnd(muB,sigma,300)';
X(:,601:900) = mvnrnd(muC,sigma,300)';

%% 3D - Cube vertices

dim_data = 100;
num_datapts = 800;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate multidimensional mean vectors
mu = zeros(dim_data,8);

mu(12:14,2) = [1,0,0];
mu(12:14,3) = [0,1,0];
mu(12:14,4) = [0,0,1];
mu(12:14,5) = [1,0,1];
mu(12:14,6) = [1,1,0];
mu(12:14,7) = [0,1,1];
mu(12:14,8) = [1,1,1];

mu = 10*mu;
% generate multidimensional covariance matrix
sigma = diag(mvnrnd(zeros(dim_data,1),eye(dim_data),1));
sigma = abs(sigma);  %enforce symmetric positive definite

for j = 1:8
    X(:,((j-1)*100 + 1):(j)*100) = mvnrnd(mu(:,j),sigma,100)';
end

%% 3 vector cube corners

dim_data = 100;
num_datapts = 900;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate multidimensional cube vectors
mu = mvnrnd(zeros(dim_data,1),eye(dim_data),3);

for j = 1:900
    X(:,j) = sign(-1 + 2*rand(1,1))*mu(1,:)' + ...
       sign(-1 + 2*rand(1,1))*mu(2,:)' + ... 
       sign(-1 + 2*rand(1,1))*mu(3,:)' + 0.15*rand(dim_data,1);
end

%% 3 vector cube

dim_data = 100;
num_datapts = 900;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate multidimensional cube vectors
mu = mvnrnd(zeros(dim_data,1),eye(dim_data),3);

% create a convex combinaion of the vectors
for k = 1:900
    lambda1 = rand(1,1);
    lambda2 = rand(1,1);
    while((lambda1 + lambda2) >1)
        lambda1 = rand(1,1);
        lambda2 = rand(1,1);
    end
    lambda3 = 1 - lambda1 - lambda2;
    
    X(:,k) = lambda1*mu(1,:)' + ...
        lambda2*mu(2,:)' + ...
        lambda3*mu(3,:)';
end



%% MNIST
addpath('DATASETS/MNIST');

% trimgid = fopen('train-images.idx3-ubyte');    
% tsimgid = fopen('t10k-images.idx3-ubyte');    
% 
% % read train data    
% fread(trimgid, 4);    
% numtrimg = toint(fread(trimgid, 4));    
% trimgh = toint(fread(trimgid, 4));    
% trimgw = toint(fread(trimgid, 4));    
% trainimages = permute(reshape(fread(trimgid,trimgh*trimgw*numtrimg),trimgh,trimgw,numtrimg), [2 1 3]);
% % read test data    
% fread(tsimgid, 4);    
% numtsimg = toint(fread(tsimgid, 4));    
% tsimgh = toint(fread(tsimgid, 4));    
% tsimgw = toint(fread(tsimgid, 4));    
% testimages = permute(reshape(fread(tsimgid, tsimgh*tsimgw*numtsimg),tsimgh,tsimgw,numtsimg), [2 1 3]);
% 


% [Test , label] = readMNIST('DATASETS/MNIST/t10k-images-idx3-ubyte','DATASETS/MNIST/t10k-labels-idx1-ubyte',1000,0);

oldpath = addpath(fullfile(matlabroot,'examples','nnet','main'));
filenameImagesTrain = 'DATASETS/MNIST/train-images.idx3-ubyte';
filenameLabelsTrain = 'DATASETS/MNIST/train-labels.idx1-ubyte';
filenameImagesTest = 'DATASETS/MNIST/t10k-images.idx3-ubyte';
filenameLabelsTest = 'DATASETS/MNIST/t10k-labels.idx1-ubyte';

XTrain = processImagesMNIST(filenameImagesTrain);
YTrain = processLabelsMNIST(filenameLabelsTrain);
XTest = processImagesMNIST(filenameImagesTest);
YTest = processLabelsMNIST(filenameLabelsTest);

%% square with corners (1,1), (1,2), (2,1), (2,2)

dim_data = 2;
num_datapts = 900;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate multidimensional cube vectors
% mu = mvnrnd(zeros(dim_data,1),eye(dim_data),3);

for k = 1:900
    X(:,k) = [0;1]*rand(1,1) + [1;0]*rand(1,1) + [1;1];
end

%% 3 orthogonal vector cube

dim_data = 100;
num_datapts = 900;

X = zeros(dim_data,num_datapts);

rng('default');  % seed random number generator

% generate orthogonal cube vectors
seed = rand(100,100);
[Q,R] = qr(seed);
mu = Q(1:3,:);

% create a convex combinaion of the vectors
for k = 1:900
    lambda1 = rand(1,1);
    lambda2 = rand(1,1);
    while((lambda1 + lambda2) >1)
        lambda1 = rand(1,1);
        lambda2 = rand(1,1);
    end
    lambda3 = 1 - lambda1 - lambda2;
    
    X(:,k) = lambda1*mu(1,:)' + ...
        lambda2*mu(2,:)' + ...
        lambda3*mu(3,:)';
end

