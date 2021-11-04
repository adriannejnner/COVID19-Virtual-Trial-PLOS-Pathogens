function mild_neut_changes_plotting(timeM,solM,timeS,solS,timeR,solR,p)

load('Human_viral_load_data.mat')

% Model curves -----------------------------------------------------------

figure
hold on 
plot(time_s5,viral_load_s5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_s6,viral_load_s6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_s18,viral_load_s18,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_g1,viral_load_g1,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_g2,viral_load_g2,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_g5,viral_load_g5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_g6,viral_load_g6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time_g7,viral_load_g7,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(timeM, solM(1,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(1,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(1,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
set(gca,'yscale','linear')
title('Viral load')
ylabel('log_{10}(copies/ml)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
legend('Data','Model')
saveas(gcf,'Fig_7A.fig');
saveas(gcf,'Fig_7A.png');


fig = figure;
hold on 
l4 = plot(timeM, (solM(2,:)+solM(4,:))*1e9,'Color',[32 52 79]/255,'LineWidth',3);
l5 = plot(timeS, (solS(2,:)+solS(4,:))*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, (solR(2,:)+solR(4,:))*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
set(gca,'yscale','log')
ylabel('cells/ml')
title('Uninfected cells')
%legend([l1 l2 l3 l4 l5], {'V(t) (mild)','Data','V(t) (severe)','S(t)+R(t) (mild)','S(t)+R(t) (severe)'})
set(gca,'Fontsize',24)
set(gca,'yTick',[1e5 1e6 1e7 1e8])
xlabel('Time (days)')
xlim([0 20])
saveas(gcf,'Fig_7B.fig');
saveas(gcf,'Fig_7B.png');
 
%dead cells
fig = figure;
hold on 
plot(timeM, solM(3,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(3,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(3,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('cells/ml')
title('Infected cells')
set(gca,'yscale','log')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_7C.fig');
saveas(gcf,'Fig_7C.png');

fig = figure
hold on 
plot(timeM, solM(5,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(5,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(5,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
set(gca,'yscale','log')
ylabel('cells/ml')
title('Dead cells')
%legend('I(t) (mild)','I(t) (severe)','D(t) (mild)','D(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1E4 1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_7D.fig');
saveas(gcf,'Fig_7D.png');

%macs cells
fig = figure;
hold on 
plot(timeM, solM(6,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(6,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(6,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('cells/ml')
title('Alveolar macs')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI7E.fig');
saveas(gcf,'Fig_SI7E.png');

figure
hold on
plot(timeM, solM(7,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(7,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(7,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
%set(gca,'yscale','log')
ylabel('cells/ml')
title('Inflammatory macs')
%legend('M_{\Phi R}(t) (mild)','M_{\Phi R}(t) (severe)','M_{\Phi I}(t) (mild)','M_{\Phi I}(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7F.fig');
saveas(gcf,'Fig_7F.png');

%monocytes cells
fig = figure;
hold on 
plot(timeM, solM(8,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(8,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(8,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('cells/ml')
title('Monocytes')
set(gca,'yscale','log')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7G.fig');
saveas(gcf,'Fig_7G.png');

figure
hold on 
plot(timeM, solM(9,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(9,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(9,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
set(gca,'yscale','linear')
ylabel('cells/ml')
title('Neutrophils')
%legend('M(t) (mild)','M(t) (severe)','N(t) (mild)','N(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
%set(gca,'yTick',[1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_7H.fig');
saveas(gcf,'Fig_7H.png');

%T cells
fig = figure;
hold on 
%yyaxis left
plot(timeM, solM(10,:)*1e9,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(10,:)*1e9,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(10,:)*1e9,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('cells/ml')
title('CD8^+ T cells')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
%legend('T(t) (mild)','T(t) (severe)')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7I.fig');
saveas(gcf,'Fig_7I.png');

%IL-6 cells
fig = figure;
hold on 
plot(timeM, solM(11,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(11,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(11,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Unbound IL-6')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7J.fig');
saveas(gcf,'Fig_7J.png');

figure
hold on 
plot(timeM, solM(12,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(12,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(12,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Bound IL-6')
%legend('L_U(t) (mild)','L_U(t) (severe)','L_B(t) (mild)','L_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI7K.fig');
saveas(gcf,'Fig_SI7K.png');

% GM-CSF cells
fig = figure;
hold on 
plot(timeM, solM(13,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(13,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(13,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Unbound GM-CSF')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7L.fig');
saveas(gcf,'Fig_7L.png');

figure
hold on
plot(timeM, solM(14,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(14,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(14,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Bound GM-CSF')
%legend('G_U(t) (mild)','G_U(t) (severe)','G_B(t) (mild)','G_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI7M.fig');
saveas(gcf,'Fig_SI7M.png');

% G-CSF cells
fig = figure;
hold on 
plot(timeM, solM(15,:)*1000,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(15,:)*1000,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(15,:)*1000,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Unbound G-CSF')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7N.fig');
saveas(gcf,'Fig_7N.png');

figure
hold on
plot(timeM, solM(16,:)*1000,'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(16,:)*1000,'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(16,:)*1000,':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Bound G-CSF')
%legend('C_U(t) (mild)','C_U(t) (severe)','C_B(t) (mild)','C_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI7O.fig');
saveas(gcf,'Fig_SI7O.png');

% IFN cells
fig = figure;
hold on 
plot(timeM, solM(17,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeS, solS(17,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
plot(timeR, solR(17,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
ylabel('pg/ml')
title('Unbound IFN ')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7P.fig');
saveas(gcf,'Fig_7P.png');

figure
hold on
plot(timeM, solM(18,:),'Color',[32 52 79]/255,'LineWidth',3);
plot(timeR, solR(18,:),':','Color',[0.47,0.64,0.19],'LineWidth',3);
plot(timeS, solS(18,:),'--','Color',[1,0.41,0.16],'LineWidth',3);
ylabel('pg/ml')
title('Bound IFN ')
%legend('F_U(t) (mild)','F_U(t) (severe)','F_B(t) (mild)','F_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')

legend('Mild disease','Neutrophils reduced','Neutrophils knockout')

saveas(gcf,'Fig_SI7Q.fig');
saveas(gcf,'Fig_SI7Q.png');

end