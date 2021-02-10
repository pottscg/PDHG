function AA_Perturbation_Trials(X, filename, noc)
    addpath('PCHA/');
    [m,n] = size(X);
    
    delta=0;
    opts.maxiter=1000;
    opts.conv_crit=1e-6;
    U = 1:n;
    I = 1:n;

    [XC,S,C,SSE,varexpl]=PCHA(X,noc,I,U,delta,opts);
    
    save(strcat('RESULTS/PERTURBATION_TRIALS/Morup_',filename,'_AA', num2str(noc)));
    
    L2_Morup_Recon = (norm(X-XC*S, 'fro'))^2*ones(1,100);
    L2_Perturbed_A_Recon = ones(1,100);
    L2_Perturbed_B_Recon = ones(1,100);
    L2_Perturbed_XB_Recon = ones(1,100);
    
    for k = 1:100
        err = k/100;
        % add random perturbation
        pert = zeros(size(S));
        pert(:,1) = err.*(-1+2.*rand(size(pert,1),1));
        A_pert = S + pert;
        
        pert = zeros(size(C));
        pert(1,:) = err.*(-1+2.*rand(1,size(pert,2)));
        B_pert = C + pert;
        
        pert = zeros(size(XC));
        pert(1,:) = err.*(-1+2.*rand(1,size(pert,2)));
        XB_pert = XC + pert;

        % Calculate reconstructions
        Recon_A_pert = X*C*A_pert;
        Recon_B_pert = X*B_pert*S;
        Recon_XB_pert = XB_pert*S;
        
        % Save result
        L2_Perturbed_A_Recon(k) = norm(X - Recon_A_pert, 'fro')^2;
        L2_Perturbed_B_Recon(k) = norm(X - Recon_B_pert, 'fro')^2;
        L2_Perturbed_XB_Recon(k) = norm(X - Recon_XB_pert, 'fro')^2;
    end
    
    figure; hold on;
    plot((1:100)./100, L2_Morup_Recon);
    plot((1:100)./100, L2_Perturbed_A_Recon, 'Color', [0.85, 0.325, 0.098]);
    hold off; 
    legend({'Morup','Perturbed A'});
    xlabel('Perturbation Strength');
    ylabel('L2 Error');
    title('MORUP VS Perturbed A Reconstruction Error');
    hold off;

    saveas(gcf,strcat('FIGURES/PERTURBATION_TRIALS/perturbed_A_plot_',filename,'_AA', num2str(noc)),'png');
    
    figure; hold on;
    plot((1:100)./100, L2_Morup_Recon);
    plot((1:100)./100, L2_Perturbed_B_Recon, 'Color', [0.85, 0.325, 0.098]);
    hold off; 
    legend({'Morup','Perturbed B'});
    xlabel('Perturbation Strength');
    ylabel('L2 Error');
    title('MORUP VS Perturbed B Reconstruction Error');
    hold off;

    saveas(gcf,strcat('FIGURES/PERTURBATION_TRIALS/perturbed_B_plot_',filename,'_AA', num2str(noc)),'png');
    
    figure; hold on;
    plot((1:100)./100, L2_Morup_Recon);
    plot((1:100)./100, L2_Perturbed_XB_Recon, 'Color', [0.85, 0.325, 0.098]);
    hold off; 
    legend({'Morup','Perturbed XB'});
    xlabel('Perturbation Strength');
    ylabel('L2 Error');
    title('MORUP VS Perturbed XB Reconstruction Error');
    hold off;

    saveas(gcf,strcat('FIGURES/PERTURBATION_TRIALS/perturbed_XB_plot_',filename,'_AA', num2str(noc)),'png');
    
    save(strcat('RESULTS/PERTURBATION_TRIALS/sketch_trials_',filename,'_AA', num2str(noc)), 'L2_Morup_Recon', 'L2_Perturbed_A_Recon', 'L2_Perturbed_B_Recon', 'L2_Perturbed_XB_Recon');
    
end
