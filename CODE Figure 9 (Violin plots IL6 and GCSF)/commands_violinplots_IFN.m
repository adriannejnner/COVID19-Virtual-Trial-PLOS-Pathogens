
load('Data_Trouillet_Assant.mat')
%epithelial cells vs alveolar macs
figure
hold on 
fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IFN Trouillet_Assant_healthy_max_IFN 1e-3],[0.9 0.9 0.9])
%plot(time_deval, solM(17,:),'Color',[32 52 79]/255,'LineWidth',3);
%plot(time_deval, solS(17,:),'--','Color',[233 150 156]/255,'LineWidth',3);
plot(data_points_Trouillet_Assant_IFN(:,1)*1000,data_points_Trouillet_Assant_IFN(:,2),'o','Color',[0.5 0.5 0.5],'LineWidth',1)
ylabel('IFN unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('Healthy range','F_U(t) (mild)','F_U(t) (severe)','IFN-postive patient')
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')
ylim([1e-3 1e2])
xlim([0 30])

for ii = 1:15
   locs =  find(data_points_Trouillet_Assant_IFN(:,1)==ii);
   mean_vec(ii) = mean(data_points_Trouillet_Assant_IFN(locs,2));
end

figure
hold on
fill([4 7 10 15 15 10 7 4],[[33,30 7 1]+[20 30 6 0],[0.002, 0.002, 0.09, 0.205]-[0 0 0.03 0.03]],[206,231,227]/255)
plot([4 7 10 15],[33,30 7 1]+[20 30 6 0],'o-','Color',[34,70,63]/255,'LineWidth',2,'MarkerSize',6)
plot([15 10 7 4],[0.002, 0.002, 0.09, 0.205]-[0 0 0.03 0.03],'o-','Color',[34,70,63]/255,'LineWidth',2,'MarkerSize',6)
scatter(data_points_Trouillet_Assant_IFN(:,1)*1000,data_points_Trouillet_Assant_IFN(:,2),30,[137 197 186]/255,'filled')
xlabel('Time (days)')
ylabel('pg/ml')
title('IFN')
set(gca,'FontSize',20)
set(gca,'yscale','log')
xlim([4 15])