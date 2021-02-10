filename = '100dim_900pt_3orthogvectorcubeconvexinterior';
load(strcat('DATASETS/',filename,'.mat'));

parfor noc = 2:50
    Sketch_Reconstruction_Trials(X,filename,noc);
end