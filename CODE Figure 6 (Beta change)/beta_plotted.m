function beta_plotted(time0,sol0,time10,sol10,time20,sol20,time30,sol30,time40,sol40,time50,sol50,p)

col = [199,233,180
127,205,187
65,182,196
29,145,192
34,94,168
12,44,132]/255;

% Model curves -----------------------------------------------------------

fig = figure;
hold on 
plot(time0,sol0(1,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(1,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(1,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(1,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(1,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(1,:)*1e9,'Color',col(1,:),'LineWidth',2)
set(gca,'yscale','linear')
title('Viral load')
ylabel('log_{10}(copies/ml)')
%legend([l1 l2 l3], {'V(t) (mild)','Data','V(t) (severe)'})
set(gca,'Fontsize',24)
xlabel('Time (days)')
xlim([0 20])
saveas(gcf,'Fig_6A.fig');
saveas(gcf,'Fig_6A.png');


fig = figure;
hold on 
plot(time0,(sol0(2,:)+sol0(4,:))*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,(sol10(2,:)+sol10(4,:))*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,(sol20(2,:)+sol20(4,:))*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,(sol30(2,:)+sol30(4,:))*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,(sol40(2,:)+sol40(4,:))*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,(sol50(2,:)+sol50(4,:))*1e9,'Color',col(1,:),'LineWidth',2)
set(gca,'yscale','log')
ylabel('cells/ml')
title('Uninfected cells')
%legend([l1 l2 l3 l4 l5], {'V(t) (mild)','Data','V(t) (severe)','S(t)+R(t) (mild)','S(t)+R(t) (severe)'})
set(gca,'Fontsize',24)
set(gca,'yTick',[1e5 1e6 1e7 1e8])
xlabel('Time (days)')
xlim([0 20])
saveas(gcf,'Fig_6B.fig');
saveas(gcf,'Fig_6B.png');
 
%dead cells
fig = figure;
hold on 
plot(time0,sol0(3,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(3,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(3,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(3,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(3,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(3,:)*1e9,'Color',col(1,:),'LineWidth',2)
ylabel('cells/ml')
title('Infected cells')
set(gca,'yscale','log')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_6C.fig');
saveas(gcf,'Fig_6C.png');

fig = figure
hold on 
plot(time0,sol0(5,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(5,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(5,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(5,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(5,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(5,:)*1e9,'Color',col(1,:),'LineWidth',2)
set(gca,'yscale','log')
ylabel('cells/ml')
title('Dead cells')
%legend('I(t) (mild)','I(t) (severe)','D(t) (mild)','D(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1E4 1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_6D.fig');
saveas(gcf,'Fig_6D.png');

%macs cells
fig = figure;
hold on 
plot(time0,sol0(6,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(6,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(6,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(6,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(6,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(6,:)*1e9,'Color',col(1,:),'LineWidth',2)
ylabel('cells/ml')
title('Alveolar macs')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI6_E.fig');
saveas(gcf,'Fig_SI6_E.png');

figure
hold on
plot(time0,sol0(7,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(7,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(7,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(7,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(7,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(7,:)*1e9,'Color',col(1,:),'LineWidth',2)
set(gca,'yscale','log')
ylabel('cells/ml')
title('Inflammatory macs')
%legend('M_{\Phi R}(t) (mild)','M_{\Phi R}(t) (severe)','M_{\Phi I}(t) (mild)','M_{\Phi I}(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6F.fig');
saveas(gcf,'Fig_6F.png');

%monocytes cells
fig = figure;
hold on 
plot(time0,sol0(8,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(8,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(8,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(8,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(8,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(8,:)*1e9,'Color',col(1,:),'LineWidth',2)
ylabel('cells/ml')
title('Monocytes')
set(gca,'yscale','log')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6G.fig');
saveas(gcf,'Fig_6G.png');

figure
hold on 
plot(time0,sol0(9,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(9,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(9,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(9,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(9,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(9,:)*1e9,'Color',col(1,:),'LineWidth',2)
set(gca,'yscale','log')
ylabel('cells/ml')
title('Neutrophils')
%legend('M(t) (mild)','M(t) (severe)','N(t) (mild)','N(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
set(gca,'yTick',[1e5 1e6 1e7 1e8])
saveas(gcf,'Fig_6H.fig');
saveas(gcf,'Fig_6H.png');

%T cells
fig = figure;
hold on 
%yyaxis left
plot(time0,sol0(10,:)*1e9,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(10,:)*1e9,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(10,:)*1e9,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(10,:)*1e9,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(10,:)*1e9,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(10,:)*1e9,'Color',col(1,:),'LineWidth',2)
ylabel('cells/ml')
title('CD8^+ T cells')
%yyaxis right
%plot(timeM,1./(1+solM(12,:)./p.eps_L_T),'Color',[230,97,1]/255,'LineWidth',2)
%plot(timeS,1./(1+solS(12,:)./p.eps_L_T),'Color',[230,97,13]/255,'LineWidth',2)
%ylabel('IL-6 anti-inflam effect, 1/(1+L_B/\epsilon_{L,T})')
%legend('T(t) (mild)','T(t) (severe)')%,'IL-6 anti-inflam curve (mild)','IL-6 anti-inflam curve (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6I.fig');
saveas(gcf,'Fig_6I.png');

%IL-6 cells
fig = figure;
hold on 
plot(time0,sol0(11,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(11,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(11,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(11,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(11,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(11,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Unbound IL-6')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6J.fig');
saveas(gcf,'Fig_6J.png');

figure
hold on 
plot(time0,sol0(12,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(12,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(12,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(12,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(12,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(12,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Bound IL-6')
%legend('L_U(t) (mild)','L_U(t) (severe)','L_B(t) (mild)','L_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI6K.fig');
saveas(gcf,'Fig_SI6K.png');

% GM-CSF cells
fig = figure;
hold on 
plot(time0,sol0(13,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(13,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(13,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(13,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(13,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(13,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Unbound GM-CSF')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6L.fig');
saveas(gcf,'Fig_6L.png');

figure
hold on
plot(time0,sol0(14,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(14,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(14,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(14,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(14,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(14,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Bound GM-CSF')
%legend('G_U(t) (mild)','G_U(t) (severe)','G_B(t) (mild)','G_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI6M.fig');
saveas(gcf,'Fig_SI6M.png');

% G-CSF cells
fig = figure;
hold on 
plot(time0,sol0(15,:)*1000,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(15,:)*1000,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(15,:)*1000,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(15,:)*1000,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(15,:)*1000,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(15,:)*1000,'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Unbound G-CSF')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6N.fig');
saveas(gcf,'Fig_6N.png');

figure
hold on
plot(time0,sol0(16,:)*1000,'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(16,:)*1000,'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(16,:)*1000,'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(16,:)*1000,'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(16,:)*1000,'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(16,:)*1000,'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Bound G-CSF')
%legend('C_U(t) (mild)','C_U(t) (severe)','C_B(t) (mild)','C_B(t) (severe)')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI6O.fig');
saveas(gcf,'Fig_SI6O.png');

% IFN cells
fig = figure;
hold on 
plot(time0,sol0(17,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(17,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(17,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(17,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(17,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(17,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Unbound IFN ')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_6P.fig');
saveas(gcf,'Fig_6P.png');

figure
hold on
plot(time0,sol0(18,:),'Color',col(6,:),'LineWidth',2)
plot(time10,sol10(18,:),'Color',col(5,:),'LineWidth',2)
plot(time20,sol20(18,:),'Color',col(4,:),'LineWidth',2)
plot(time30,sol30(18,:),'Color',col(3,:),'LineWidth',2)
plot(time40,sol40(18,:),'Color',col(2,:),'LineWidth',2)
plot(time50,sol50(18,:),'Color',col(1,:),'LineWidth',2)
ylabel('pg/ml')
title('Bound IFN ')
%legend('F_U(t) (mild)','F_U(t) (severe)','F_B(t) (mild)','F_B(t) (severe)')
legend('Mild disease','Severe disease')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_SI6Q.fig');
saveas(gcf,'Fig_SI6Q.png');

figure
hold on
plot(time0,sol0(18,:),'Color',col(6,:),'LineWidth',4)
plot(time10,sol10(18,:),'Color',col(5,:),'LineWidth',4)
plot(time20,sol20(18,:),'Color',col(4,:),'LineWidth',4)
plot(time30,sol30(18,:),'Color',col(3,:),'LineWidth',4)
plot(time40,sol40(18,:),'Color',col(2,:),'LineWidth',4)
plot(time50,sol50(18,:),'Color',col(1,:),'LineWidth',4)
legend('0%','10%','20%','30%','40%','50%')
set(gca,'FontSize',24)


end