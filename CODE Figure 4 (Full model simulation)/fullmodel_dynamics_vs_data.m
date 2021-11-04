function fullmodel_dynamics_vs_data(timeM,solM,timeS,solS,p)

load('Data_Trouillet_Assant.mat')
load('Laing_data.mat')
% 
% figure
% hold on 
% fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IFN Trouillet_Assant_healthy_max_IFN 1e-3],[0.85 0.85 0.85])
% plot(timeM, solM(17,:),'Color',[94,60,153]/255,'LineWidth',2)
% plot(timeS, solS(17,:),'--','Color',[94,60,153]/255,'LineWidth',2)
% plot(data_points_Trouillet_Assant_IFN(:,1)*1000,data_points_Trouillet_Assant_IFN(:,2),'o','Color',[0.5 0.5 0.5],'LineWidth',1)
% ylabel('IFN unbound (pg/ml)')
% hold on
% set(gca,'Fontsize',18)
% xlabel('Days')
% legend('Healthy range','IFN unbound (mild), F_U(t)','IFN unbound (severe), F_U(t)','IFN-postive SARS-CoV-2 infection')
% ax = gca;
% ax.YAxis(1).Color = [94,60,153]/255;
% set(gca,'yscale','log')
% ylim([1e-3 1e2])
% 
% figure
% hold on 
% fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IL6 Trouillet_Assant_healthy_max_IL6 1e-3],[0.85 0.85 0.85])
% plot(timeM, solM(11,:),'Color',[94,60,153]/255,'LineWidth',2)
% plot(timeS, solS(11,:),'--','Color',[94,60,153]/255,'LineWidth',2)
% plot(data_points_Trouillet_Assant_IL6_IFNpos(:,1),data_points_Trouillet_Assant_IL6_IFNpos(:,2),'o','Color',[0.5 0.5 0.5],'LineWidth',1)
% plot(data_points_Trouillet_Assant_IL6_IFNneg(:,1),data_points_Trouillet_Assant_IL6_IFNneg(:,2),'*','Color',[0.5 0.5 0.5],'LineWidth',1)
% ylabel('IL-6 unbound (pg/ml)')
% hold on
% set(gca,'Fontsize',18)
% xlabel('Days')
% legend('Healthy range','L_U(t) (mild)','L_U(t) (severe)','IFN-postive','IFN-negative')
% ax = gca;
% ax.YAxis(1).Color = [94,60,153]/255;
% set(gca,'yscale','log')
% ylim([1e-1 1e5])

load('IL6_Herold.mat')
figure
hold on
fill([0 0 30 30],[IL6_ranges_mild(1) IL6_ranges_mild(2) IL6_ranges_mild(2) IL6_ranges_mild(1)],[0.85 0.85 0.85])
fill([0 0 30 30],[IL6_ranges_severe(1) IL6_ranges_severe(2) IL6_ranges_severe(2) IL6_ranges_severe(1)],[0.85 0.85 0.85])
plot(timeM, solM(11,:),'Color',[94,60,153]/255,'LineWidth',2)
plot(timeS, solS(11,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('IL-6 unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('No mechanical ventillation','Yes mechanical ventillation','L_U(t) (mild)','L_U(t) (severe)')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')

load('Long_measurements.mat')
figure
hold on
fill([0 0 30 30],[1E-1 IL6_Long(1) IL6_Long(1) 1E-1],[0.85 0.85 0.85])
fill([0 0 30 30],[1E-1 IL6_Long(2) IL6_Long(2) 1E-1],[0.85 0.85 0.85])
plot(timeM, solM(11,:),'Color',[94,60,153]/255,'LineWidth',2)
plot(timeS, solS(11,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('IL-6 unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('AS','S','L_U(t) (mild)','L_U(t) (severe)')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','linear')

figure
hold on
fill([0 0 30 30],[1e-6 GCSF_Long(1)/1000 GCSF_Long(1)/1000 1e-6],[0.85 0.85 0.85])
fill([0 0 30 30],[1e-6 GCSF_Long(2)/1000 GCSF_Long(2)/1000 1e-6],[0.85 0.85 0.85])
plot(timeM, solM(15,:),'Color',[94,60,153]/255,'LineWidth',2)
plot(timeS, solS(15,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('G-CSF unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('AS','S','C_U(t) (mild)','C_U(t) (severe)')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')



load('Laing_data.mat')

figure
hold on 
%T cells
fig = figure;
hold on 
%yyaxis left
plot(timeM, solM(10,:)*1e9,'Color',[94,60,153]/255,'LineWidth',2);
plot(timeS, solS(10,:)*1e9,'--','Color',[94,60,153]/255,'LineWidth',2);
ylabel('CD8^+ T cells (cells/mL)')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
set(gca,'Fontsize',19)
xlabel('Time (days)')
plot(CD8_T_cells_Laing(1:10,1)+32,CD8_T_cells_Laing(1:10,2),'o','LineWidth',2)
plot(CD8_T_cells_Laing(11:11+24,1)+34,CD8_T_cells_Laing(11:11+24,2),'o','LineWidth',2)
plot(CD8_T_cells_Laing(11+24:end,1)+35,CD8_T_cells_Laing(11+24:end,2),'o','LineWidth',2)
set(gca,'yscale','log')
legend('T(t) (mild)','T(t) (severe)','Low','Moderate','Severe')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'XTick',[0 10 20 30])



figure
hold on 
%T cells
fig = figure;
hold on 
%yyaxis left
plot(timeM, solM(9,:)*1e9,'Color',[230,97,1]/255,'LineWidth',2);
plot(timeS, solS(9,:)*1e9,'--','Color',[230,97,1]/255,'LineWidth',2);
set(gca,'yscale','log')
ylabel('Neutrophils (cells/mL)')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
set(gca,'Fontsize',19)
xlabel('Time (days)')
plot(Neutrophils_Laing(1:9,1)+32,Neutrophils_Laing(1:9,2),'o')
plot(Neutrophils_Laing(10:10+19,1)+34,Neutrophils_Laing(10:10+19,2),'o')
plot(Neutrophils_Laing(10+19:end,1)+36,Neutrophils_Laing(10+19:end,2),'o')
set(gca,'yscale','log')
legend('T(t) (mild)','T(t) (severe)','Low','Moderate','Severe')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'XTick',[0 10 20 30])


figure
hold on 
%T cells
fig = figure;
hold on 
%yyaxis left
intSolM = interp1(timeM,solM(9,:),linspace(0,timeM(end),1000));
plot(linspace(0,timeM(end),1000), intSolM*1e9/min(solM(9,:)*1e9),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(9,:)*1e9/min(solM(9,:)*1e9),'--','Color',[233 150 156]/255,'LineWidth',3);
 set(gca,'yscale','linear')
ylabel('Neutrophils (cells/mL)')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')

legend('T(t) (mild)','T(t) (severe)','Low','Moderate','Severe')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'XTick',[0 10 20 30])
set(gca,'Fontsize',19)
xlabel('Time (days)')

ylim([0 9])

data1 = [Neutrophils_Laing(1:9,2)/min(Neutrophils_Laing(1:9,2))];
data2 = [Neutrophils_Laing(10:10+19,2)/min(Neutrophils_Laing(1:9,2))];
data3 = [Neutrophils_Laing(10+19:end,2)/min(Neutrophils_Laing(1:9,2))];

group = [30*ones(size(data1));31*ones(size(data2)); 32*ones(size(data3))];

figure
boxplot([data1;data2;data3],group)
set(gca,'Fontsize',19)

ylim([0 9])



figure
hold on 
%T cells
fig = figure;
hold on 
%yyaxis left
plot(timeM, solM(10,:)*1e9,'Color',[32 52 79]/255,'LineWidth',2);
plot(timeS, solS(10,:)*1e9,'--','Color',[233 150 156]/255,'LineWidth',2);
ylabel('CD8^+ T cells (cells/mL)')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
set(gca,'Fontsize',19)
xlabel('Time (days)')
%legend('T(t) (mild)','T(t) (severe)','Low','Moderate','Severe')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'XTick',[0 10 20 30])


data1 = [CD8_T_cells_Laing(1:10,2)];
data2 = [CD8_T_cells_Laing(11:11+24,2)];
data3 = [CD8_T_cells_Laing(11+24:end,2)];

group = [30*ones(size(data1));31*ones(size(data2)); 32*ones(size(data3))];

figure
boxplot([data1;data2;data3],group)
set(gca,'Fontsize',19)

end

