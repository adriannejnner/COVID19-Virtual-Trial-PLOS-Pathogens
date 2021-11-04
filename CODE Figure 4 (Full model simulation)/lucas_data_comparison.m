function lucas_data_comparison(timeM,solM,timeS,solS,p)

load('Lucas_data.mat')
%monocytes cells
fig = figure;
hold on 
plot(timeM, solM(8,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(8,:)*1e9,'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('cells/ml')
title('Monocytes')
set(gca,'yscale','log')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_4G_2.fig');
saveas(gcf,'Fig_4G_2.png');

violin_mat_monocyte = [[Monocytes_Lucas_HCW;NaN(48-28,1)],Monocytes_Lucas_moderate,[Monocytes_Lucas_Severe;NaN(48-29,1)]];
figure
hold on 
violinplot(violin_mat_monocyte)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('Monocytes (cells/ml)')

figure
hold on 
plot(timeM, solM(9,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(9,:)*1e9,'--','Color',[233 150 156]/255,'LineWidth',3);
set(gca,'yscale','log')
ylabel('cells/ml')
title('Neutrophils')
%legend('M(t) (mild)','M(t) (severe)','N(t) (mild)','N(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_4H_2.fig');
saveas(gcf,'Fig_4H_2.png');


violin_mat_neutrophil = [[Neutrophils_Lucas_HCW;NaN(57-20,1)],Neutrophils_Lucas_moderate,[Neutrophils_Lucas_Severe;NaN(57-32,1)]];
figure
hold on 
violinplot(violin_mat_neutrophil)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('Neutrophils (cells/ml)')

%T cells
fig = figure;
hold on 
%yyaxis left
plot(timeM, solM(10,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(10,:)*1e9,'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('cells/ml')
title('CD8^+ T cells')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
%legend('T(t) (mild)','T(t) (severe)')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_4I_2.fig');
saveas(gcf,'Fig_4I_2.png');

violin_mat_Tcell = [[Tcells_Lucas_HCW;NaN(69-43,1)],Tcells_Lucas_moderate,[Tcells_Lucas_Severe;NaN(69-32,1)]];
figure
hold on 
violinplot(violin_mat_Tcell)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('T cells (cells/ml)')

%IL-6 cells
fig = figure;
hold on 
plot(timeM, solM(11,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(11,:),'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('pg/ml')
title('Unbound IL-6')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_4J_2.fig');
saveas(gcf,'Fig_4J_2.png');


violin_mat_IL6 = [[IL6_Lucas_HCW_individual';NaN(108-47,1)],IL6_Lucas_moderate_individual',[IL6_Lucas_Severe_individual';NaN(108-41,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('IL-6 (log(pg/mL))')



% IFN cells
fig = figure;
hold on 
plot(timeM, solM(17,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(17,:),'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('pg/ml')
title('Unbound IFN ')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_4P_2.fig');
saveas(gcf,'Fig_4P_2.png');

figure
hold on
plot(timeM, solM(18,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(18,:),'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('pg/ml')
title('Bound IFN ')
%legend('F_U(t) (mild)','F_U(t) (severe)','F_B(t) (mild)','F_B(t) (severe)')
legend('Mild disease','Severe disease')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_4Q_2.fig');
saveas(gcf,'Fig_4Q_2.png');


violin_mat_IFN = [[IFN_Lucas_HCW;NaN(102-46,1)],IFN_Lucas_moderate,[IFN_Lucas_Severe;NaN(102-46,1)]];
figure
hold on 
violinplot(violin_mat_IFN)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('IFN (pg/ml)')



end