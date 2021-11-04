%simulating full model

format long

p = load_parameters();
p.F_U_0 = 0.015;%5;
p = Homeostasis_calculations_IFN(p);

p.eps_F_I = 2*1e-4;%10*1e-5;%5;%133%650*p.F_B_0/p.F_U_0;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1*2;

p.phat =394;%741.2;
p.d_V = 8.4;%18.94;
p.beta =0.3;%0.289;
p.d_I = 0.1;%0.144;
p.V0 = 4.5;%12;

p.S0 = 0.16;
p.Smax = 0.16;



%fixed parameters
p.r = 0.7397;
p.d_D = 8;
p.Smax = 0.16;
p.S0 = 0.16;
p.I0 = 0;
p.D0 = 0;
p.tau_I = 0.1667;
p.lag_s5 = 2;
p.lag_s6 = 2.79;
p.lag_s18 = 1.32;
p.lag_g1 = 1.27;
p.lag_g2 = 0.92;
p.lag_g5 = 1.32;
p.lag_g6 = 2.76;
p.lag_g7 = 2;
p.phat = mean(10.^[2.59,2.6,2.59,2.59,2.6, 2.59, 2.6, 2.59, 2.59, 2.59, 2.59, 2.59, 2.6, 2.6 2.6 2.6 2.6 2.6 2.6 2.59, 2.59, 2.6, 2.6, 2.6]);

tspan = [0 30];

[timeO,solO] = COVID_IMMUNE_MODEL_virus_IFN_resistance(p,tspan);
%[timeO,solO] = COVID_IMMUNE_MODEL_virus_IFN_resistance_noRprolif(p,tspan);

%p.p_F_I = 0.2;
%p.eta_F_I = 0.1;
[timeD,solD] = COVID_IMMUNE_MODEL_virus_IFN_resistance_delay(p,tspan); %IFN delay
%[timeD,solD] = COVID_IMMUNE_MODEL_virus_IFN_resistance_delay_noRprolif(p,tspan); %IFN delay

figure
hold on 
yyaxis left
plot(1./(1+solO(7,:)./(p.eps_F_I)),'Color',[94,60,153]/255,'LineWidth',2)
ylabel('Effect curve 1./(1+F_B/\epsilon_{F,I})')
yyaxis right
plot(1./(1+solD(7,:)./(p.eps_F_I)),'Color',[230,97,1]/255,'LineWidth',2)
ylabel('Effect curve 1./(1+F_B/\epsilon_{F,I})')
set(gca,'Fontsize',18)
xlabel('Days')
legend('Effect curve','Effect curve-IFN delay')

ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
ax.YAxis(2).Color = [230,97,1]/255;

load('Human_viral_load_data.mat')

% Model curves -----------------------------------------------------------

figure
hold on  
l3 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
l1 = plot(timeO, solO(1,:),'Color',[94,60,153]/255,'LineWidth',2)
l2 = plot(timeD, solD(1,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('Viral load (log_{10}(copies/mL))')
yyaxis right
l4 = plot(timeO, (solO(2,:)+solO(4,:))*1e9,'Color',[230,97,1]/255,'LineWidth',2)
l5 = plot(timeD, (solD(2,:)+solD(4,:))*1e9,'--','Color',[230,97,1]/255,'LineWidth',2)
set(gca,'yscale','log')
ylabel('Cells/ml')
legend([l1 l2 l3 l4 l5], {'V(t)','V(t)-IFN delay','Viral load data','S(t)+R(t)','S(t)+R(t)-IFN delay'})
set(gca,'Fontsize',18)
xlabel('Time (days)')
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;%[230,97,1]/255;
ax.YAxis(2).Color = [230,97,1]/255;
ylim([1e6 1e10])

%plot(tspan, [2 2],'k--','LineWidth',2)
set(gca,'FontSize',18)
%legend([l1 l2],'Viral load','V(t)')
xlim(tspan)

colmap = [140,81,10;...
191,129,45;...
223,194,125;...
246,232,195;...
199,234,229;...
128,205,193
53,151,143;...
1,102,94]/255;

figure
hold on  
l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[140,81,10]/255,'LineWidth',2.5)
l2 = plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
l3 = plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
l4 = plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
l5 = plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
l6 = plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
l7 = plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
l8 = plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)
l9 = plot(timeO, solO(1,:),'Color',[94,60,153]/255,'LineWidth',2)
l10 = plot(timeD, solD(1,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('Viral load (log_{10}(copies/mL))')
yyaxis right
l11 = plot(timeO, (solO(2,:)+solO(4,:))*1e9,'Color',[230,97,1]/255,'LineWidth',2)
l12 = plot(timeD, (solD(2,:)+solD(4,:))*1e9,'--','Color',[230,97,1]/255,'LineWidth',2)
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12],{'S5','S6','S18','G1','G2','G5','G6','G7','V(t)','V(t)-IFN delay','Viral load data','S(t)+R(t)','S(t)+R(t)-IFN delay'})
plot(tspan, [2 2],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'FontSize',18)
set(gca,'yscale','log')
xlim(tspan)
xlabel('Time (days)')
ylabel('Cells/ml')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;%[230,97,1]/255;
ax.YAxis(2).Color = [230,97,1]/255;
ylim([1e6 1e10])

figure
hold on 
yyaxis left
plot(timeO, solO(3,:)*1e9,'Color',[94,60,153]/255,'LineWidth',2)
plot(timeD, solD(3,:)*1e9,'--','Color',[94,60,153]/255,'LineWidth',2)
%plot(timeO, solO(4,:)*1e9,'k-','LineWidth',2)
%plot(timeD, solD(4,:)*1e9,'k:','LineWidth',2)
set(gca,'yscale','log')
ylabel('Infected cells (cells/mL)')
yyaxis right
hold on
plot(timeO, solO(5,:)*1e9,'Color',[230,97,1]/255,'LineWidth',2)
plot(timeD, solD(5,:)*1e9,'--','Color',[230,97,1]/255,'LineWidth',2)
set(gca,'yscale','log')
ylabel('Dead cells (cells/ml)')
set(gca,'Fontsize',18)
xlabel('Time (days)')
legend('I(t)','I(t)-IFN delay','D(t)','D(t)-IFN delay')%,'R(t)','R(t)-IFN delay'
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
ax.YAxis(2).Color = [230,97,1]/255;

%epithelial cells vs alveolar macs
figure
hold on 
yyaxis left
plot(timeO, solO(6,:),'Color',[94,60,153]/255,'LineWidth',2)
plot(timeD, solD(6,:),'--','Color',[94,60,153]/255,'LineWidth',2)
ylabel('IFN unbound (pg/ml)')
yyaxis right
hold on
plot(timeO, solO(7,:),'Color',[230,97,1]/255,'LineWidth',2)
plot(timeD, solD(7,:),'--','Color',[230,97,1]/255,'LineWidth',2)
ylabel('IFN bound (pg/ml)')
set(gca,'Fontsize',18)
xlabel('Days')
legend('F_U(t)','F_U(t)-IFN delay','F_B(t)','F_B(t)-IFN delay')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
ax.YAxis(2).Color = [230,97,1]/255;

load('Data_Trouillet_Assant.mat')
%epithelial cells vs alveolar macs
figure
hold on 
fill([0 0 40 40],[1e-3 Trouillet_Assant_healthy_max_IFN Trouillet_Assant_healthy_max_IFN 1e-3],[0.9 0.9 0.9])
plot(timeO, solO(6,:),'Color',[94,60,153]/255,'LineWidth',2)
plot(data_points_Trouillet_Assant_IFN(:,1)*1000,data_points_Trouillet_Assant_IFN(:,2),'o','Color',[0.5 0.5 0.5],'LineWidth',1)
ylabel('IFN unbound (pg/ml)')
hold on
set(gca,'Fontsize',18)
xlabel('Days')
legend('Healthy range','Model simulation, F_U(t)','IFN-postive SARS-CoV-2 infection')
ax = gca;
%ax.YAxis(1).Color = [94,60,153]/255;
set(gca,'yscale','log')
ylim([1e-3 1e2])
xlim([0 30])


