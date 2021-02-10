figure; 
hold on;
for j = 1:10
    plot(data_dim, Morup_err(:,j), 'color', 'b');
    plot(data_dim, PDHG_err(:,j), 'color', [0.85, 0.325, 0.098]);
end
hold off;


%%
figure; 
hold on;
for j=1:10
    scatter(data_dim, Morup_err(:,j),2,'filled', 'b'); 
    scatter(data_dim, PDHG_err(:,j),2,[0.85, 0.325, 0.098],'filled');
end
Morup_avg = sum(Morup_err,2)./10;
PDHG_avg = sum(PDHG_err,2)./10;
plot(data_dim, Morup_avg, 'color', 'b');
plot(data_dim, PDHG_avg, 'color', [0.85, 0.325, 0.098]);

hold off; 

%%
figure; 
subplot(2,1,1);
hold on;
for j=1:10
    scatter(data_dim, Morup_err(:,j),2,'g', 'filled'); 
end
for j = 1:100
    for k = 1:10
        if(Morup_constraints(j,k) == 0)
            scatter(data_dim(j), Morup_err(j,k), 2, 'r', 'filled');
        end
    end
end
Morup_avg = sum(Morup_err,2)./10;
plot(data_dim, Morup_avg, 'color', 'b');
hold off; 
title('Morup Reconstruction Error');
ylim([100,700]);
xlabel('Data Dimension');
ylabel('Model Error');

subplot(2,1,2);
hold on; 
for j=1:10
    scatter(data_dim, PDHG_err(:,j),2,'g','filled');
end
for j = 1:100
    for k = 1:10
        if(PDHG_constraints(j,k) == 0)
            scatter(data_dim(j), PDHG_err(j,k), 2, 'r', 'filled');
        end
    end
end
PDHG_avg = sum(PDHG_err,2)./10;
plot(data_dim, PDHG_avg, 'color', 'b');
title('PDHG Reconstruction Error');
ylim([100,700]);
xlabel('Data Dimension');
ylabel('Model Error');
hold off; 