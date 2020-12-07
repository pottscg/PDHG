addpath('PCHA/');

noc=20; % Number of archetypes

U=1:size(X,2); % Entries in X used that is modelled by the AA model
I=1:size(X,2); % Entries in X used to define archetypes
% if two expensive to useall entries for I find N relevant observations by
% the following procedure:
% N=100;
% I=FurthestSum(X,N,ceil(rand*size(X,2)));

delta=0;
opts.maxiter=1000;
opts.conv_crit=1e-6;

% Use PCHA.m
[XC,S,C,SSE,varexpl]=PCHA(X,noc,I,U,delta,opts);

%% run_PCHA_PDHG

addpath('PDHG/');

noc=20; % Number of archetypes

U=1:size(X,2); % Entries in X used that is modelled by the AA model
I=1:size(X,2); % Entries in X used to define archetypes
% if two expensive to useall entries for I find N relevant observations by
% the following procedure:
% N=100;
% I=FurthestSum(X,N,ceil(rand*size(X,2)));

delta=0;
opts.maxiter=1000;
opts.conv_crit=1e-6;

% Use PCHA base with PDHG solver
[XB,A,B,SSE,varexpl]=PCHA_PDHG(X,noc,I,U,delta,opts);
