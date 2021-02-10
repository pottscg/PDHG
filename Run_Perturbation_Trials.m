filename = '100dim_900pt_3orthogvectorcubeconvexinterior';
load(strcat('DATASETS/',filename,'.mat'));

parfor noc = 2:50
    AA_Perturbation_Trials(X,filename,noc);
end