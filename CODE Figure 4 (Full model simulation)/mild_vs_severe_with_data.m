function mild_vs_severe_with_data(time_deval,solM,time,solS,p)

load('Data_Trouillet_Assant.mat')
%epithelial cells vs alveolar macs
figure
hold on 
fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IFN Trouillet_Assant_healthy_max_IFN 1e-3],[0.9 0.9 0.9])
plot(time_deval, solM(17,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(time_deval, solS(17,:),'--','Color',[233 150 156]/255,'LineWidth',3);
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

figure
hold on 
fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IL6 Trouillet_Assant_healthy_max_IL6 1e-3],[0.9 0.9 0.9])
plot(time, solM(11,:),'Color',[32 52 79]/255,'LineWidth',2)
plot(time, solS(11,:),'--','Color',[233 150 156]/255,'LineWidth',2)
plot(data_points_Trouillet_Assant_IL6_IFNpos(:,1),data_points_Trouillet_Assant_IL6_IFNpos(:,2),'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(data_points_Trouillet_Assant_IL6_IFNneg(:,1),data_points_Trouillet_Assant_IL6_IFNneg(:,2),'*','Color',[0.93 0.69 0.13],'LineWidth',2)
ylabel('IL-6 unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('Healthy range','L_U(t) (mild)','L_U(t) (severe)','IFN-postive','IFN-negative')
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')
ylim([0.4 1e5])
xlim([0 30])
set(gca,'YTick',[1e0, 1e3, 1e5])


load('IL6_herold_individual.mat')
figure
hold on
violin_mat_IL6 = [no_ventilation,[with_ventilation;repmat(NaN,10,1)]];
l1 = violinplot(violin_mat_IL6)
hold on 
l2 = plot((time/4)+3, solM(11,:),'Color',[32 52 79]/255,'LineWidth',2)
l3 = plot((time/4)+3, solS(11,:),'--','Color',[233 150 156]/255,'LineWidth',2)
ylabel('IL-6 unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend([l1(1).ScatterPlot, l1(2).ScatterPlot l2 l3],{'No mechanical ventillation','Yes mechanical ventillation','L_U(t) (mild)','L_U(t) (severe)'})
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')
set(gca,'xtick',[1 2 linspace(3,(time(end)/4)+3,4)],'xticklabels',{'No','Yes','0','10','20','30'})

set(gca,'YTick',[1e0, 1e3, 1e5])
xlabel('Days')
xlim([0 30/4+3])
ylim([0.4 1e5])



load('Lucas_data.mat')
violin_mat_IL6 = [[IL6_Lucas_HCW_individual';NaN(108-47,1)],IL6_Lucas_moderate_individual',[IL6_Lucas_Severe_individual';NaN(108-41,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
set(gca,'FontSize',18)
xlabel('Response')
set(gca,'XTick',[1,2,3],'XTickLabels',{'HCW','Moderate','Severe'})
xlim([0.5 3.5])
ylabel('IL-6 (log(pg/mL))')

load('GCSF_Long_individual.mat')
figure
hold on
violin_mat_GCSF = [[AS;NaN],S];
l1 = violinplot(violin_mat_GCSF)
hold on 
l2 = plot((time/4)+3, solM(15,:)*1000,'Color',[32 52 79]/255,'LineWidth',3);
l3 = plot((time/4)+3, solS(15,:)*1000,'--','Color',[233 150 156]/255,'LineWidth',3);
ylabel('G-CSF unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend([l1(1).ScatterPlot, l1(2).ScatterPlot l2 l3],{'Asymptomatic','Symptomatic','G_U(t) (mild)','G_U(t) (severe)'})
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','linear')
set(gca,'xtick',[1 2 linspace(3,(time(end)/4)+3,4)],'xticklabels',{'AS','S','0','10','20','30'})

set(gca,'YTick',[1e1, 1e2, 2*1e2, 3e2])
xlabel('Days')
xlim([0 30/4+3])
ylim([3 350])

load('Lucas_data.mat')
violin_mat_neutrophil = [Neutrophils_Lucas_moderate,[Neutrophils_Lucas_Severe;NaN(57-32,1)]]./mean(Neutrophils_Lucas_HCW);
figure
hold on 
l1 = violinplot(violin_mat_neutrophil)
hold on 
l2 = plot((time/4)+3, solM(9,:)/p.N0,'Color',[32 52 79]/255,'LineWidth',3);
l3 = plot((time/4)+3, solS(9,:)/p.N0,'--','Color',[233 150 156]/255,'LineWidth',3);
set(gca,'FontSize',18)
xlabel('Days')
set(gca,'XTick',[1,2 linspace(3,(time(end)/4)+3,4)],'XTickLabels',{'M','S','0','10','20','30'})
legend([l1(1).ScatterPlot, l1(2).ScatterPlot l2 l3],{'Moderate patient','Severe patient','N(t) (mild)','N(t) (severe)'})
xlim([0 30/4+3])
ylabel('Change from HCW baseline')

violin_mat_monocyte = [Monocytes_Lucas_moderate,[Monocytes_Lucas_Severe;NaN(48-29,1)]]./mean(Monocytes_Lucas_HCW);
figure
hold on 
l1 = violinplot(violin_mat_monocyte)
hold on 
l2 = plot((time/4)+3, solM(8,:)/p.M0,'Color',[32 52 79]/255,'LineWidth',3);
l3 = plot((time/4)+3, solS(8,:)/p.M0,'--','Color',[233 150 156]/255,'LineWidth',3);
set(gca,'FontSize',18)
xlabel('Days')
set(gca,'XTick',[1,2 linspace(3,(time(end)/4)+3,4)],'XTickLabels',{'M','S','0','10','20','30'})
legend([l1(1).ScatterPlot, l1(2).ScatterPlot l2 l3],{'Moderate patient','Severe patient','M(t) (mild)','M(t) (severe)'})
xlim([0 30/4+3])
ylabel('Change from HCW baseline')

violin_mat_Tcells = [Tcells_Lucas_moderate,[Tcells_Lucas_Severe;NaN(69-32,1)]]./mean(Tcells_Lucas_HCW);
figure
hold on 
l1 = violinplot(violin_mat_Tcells)
hold on 
l2 = plot((time/4)+3, solM(10,:)/p.T0,'Color',[32 52 79]/255,'LineWidth',3);
l3 = plot((time/4)+3, solS(10,:)/p.T0,'--','Color',[233 150 156]/255,'LineWidth',3);
set(gca,'FontSize',18)
xlabel('Days')
set(gca,'XTick',[1,2 linspace(3,(time(end)/4)+3,4)],'XTickLabels',{'M','S','0','10','20','30'})
legend([l1(1).ScatterPlot, l1(2).ScatterPlot l2 l3],{'Moderate patient','Severe patient','T(t) (mild)','T(t) (severe)'})
xlim([0 30/4+3])
ylabel('Change from HCW baseline')
set(gca,'ytick',[1e-1 1e0 1e1])

end
