function Sketch_Reconstruction_Trials(X, filename, noc)
    addpath('PCHA/');
    [m,n] = size(X);
    
    delta=0;
    opts.maxiter=1000;
    opts.conv_crit=1e-6;
    U = 1:n;
    I = 1:n;

    [XC,S,C,SSE,varexpl]=PCHA(X,noc,I,U,delta,opts);
    
    save(strcat('RESULTS/SKETCH_TRIALS/Morup_',filename,'_AA', num2str(noc)));
    
    L2_Morup_Recon = (norm(X-XC*S, 'fro'))^2*ones(1,99);
    L2_Sketch_Recon = ones(1,99);
    
    for k = 2:100
        % create sketch data Y
        Psi = randn(k,m);
        Y = Psi*X;
        
        U=1:size(Y,2); % Entries in Y used that is modelled by the AA model
        I=1:size(Y,2); % Entries in Y used to define archetypes
        
        % run PCHA
        [YC,Stilda,Ctilda,SSE,varexpl]=PCHA(Y,noc,I,U,delta,opts);

        % Calculate archetypes from full data
        XC_sketches = X*Ctilda;
        Recon_sketches = XC_sketches*Stilda;
        
        % Save result
        L2_Sketch_Recon(k-1) = norm(X - Recon_sketches, 'fro')^2;
    end
    
    figure; hold on;
    plot(2:100, L2_Morup_Recon);
    plot(2:100, L2_Sketch_Recon, 'Color', [0.85, 0.325, 0.098]);
    hold off; 
    legend({'Morup','Sketch'});
    xlabel('Sketch Projection Dimension');
    ylabel('L2 Error');
    title('MORUP VS Sketch Reconstruction Error');

    saveas(gcf,strcat('FIGURES/SKETCH_TRIALS/sketch_trials_plot_',filename,'_AA', num2str(noc)),'png');
    
    save(strcat('RESULTS/SKETCH_TRIALS/sketch_trials_',filename,'_AA', num2str(noc)), 'L2_Morup_Recon', 'L2_Sketch_Recon');
    
end
