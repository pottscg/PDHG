
%%
figure; 
data_dim = 3:3:3*85;
% subplot(2,1,1);
hold on;
for j=1:20
    scatter(data_dim, L2_Model_err(:,j),2,[0.85, 0.325, 0.098], 'filled'); 
end
% for j = 1:100
%     for k = 1:10
%         if(Morup_C_constraints(j,k) == 1)
%             scatter(data_dim(j), Morup_err(j,k), 2, 'r', 'filled');
%         end
%     end
% end
Model_avg = sum(L2_Model_err,2)./20;
plot(data_dim, Model_avg, 'color', 'b');
plot(data_dim, norm(X-Recon_Morup, 'fro')^2*ones(size(data_dim)), 'g');
hold off; 
% title('Sketch Reconstruction Error');
ylim([4*10^(5), 2*10^6]);
xlim([0, 256]);
set(gca, 'YScale', 'log');
xlabel('Projection Dimension');
ylabel('Model Error');

% subplot(2,1,2);
% hold on; 
% for j=1:10
%     scatter(data_dim, PDHG_err(:,j),2,'g','filled');
% end
% for j = 1:100
%     for k = 1:10
%         if(PDHG_B_constraints(j,k) == 1)
%             scatter(data_dim(j), PDHG_err(j,k), 2, 'r', 'filled');
%         end
%     end
% end
% PDHG_avg = sum(PDHG_err,2)./10;
% plot(data_dim, PDHG_avg, 'color', 'b');
% title('PDHG Reconstruction Error');
% % ylim([100,700]);
% xlabel('Data Dimension');
% ylabel('Model Error');
% hold off; 