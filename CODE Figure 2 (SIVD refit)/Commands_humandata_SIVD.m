%% fit human data fixing some parameters to their value in goyal and fiting the remaining

load('Human_viral_load_data.mat')

%set up data
data.time_s5 = time_s5;
data.time_s6 = time_s6;
data.time_s18 = time_s18;
data.time_g1 = time_g1;
data.time_g2 = time_g2;
data.time_g5 = time_g5;
data.time_g6 = time_g6;
data.time_g7 = time_g7;

data.data_s5 = viral_load_s5;
data.data_s6 = viral_load_s6;
data.data_s18 = viral_load_s18;
data.data_g1 = viral_load_g1;
data.data_g2 = viral_load_g2;
data.data_g5 = viral_load_g5;
data.data_g6 = viral_load_g6;
data.data_g7 = viral_load_g7;


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
p.p = mean(10.^[2.59,2.6,2.59,2.59,2.6, 2.59, 2.6, 2.59, 2.59, 2.59, 2.59, 2.59, 2.6, 2.6 2.6 2.6 2.6 2.6 2.6 2.59, 2.59, 2.6, 2.6, 2.6])
p.max_time = max([time_s5+p.lag_s5,data.time_s6+p.lag_s6,time_s18+p.lag_s18,time_g1+p.lag_g1,time_g2+p.lag_g2,time_g5+p.lag_g5,time_g6+p.lag_g6,time_g7+p.lag_g7])

%initial parameter guesses
%p.p = 741;
p.d_V = 18.81;
p.beta = 0.3;
p.d_I = 0.14;
p.V0 = 1;

[t_av, sol_av, sol_S, sol_I, sol_D, param_fit] = fit_SIVD_parameters(p,data);

betafit = param_fit(1);
d_Vfit = param_fit(2);
d_Ifit = param_fit(3);
V0fit = param_fit(4);


figure
hold on  
l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
l2 = plot(t_av,sol_av,'Color',[230,97,1]/255,'LineWidth',2)
plot([0 p.max_time], [2 2],'k--','LineWidth',2)
set(gca,'xtick',[0 p.max_time/2,p.max_time],'xticklabel',{'S','S+10','S+20'})
set(gca,'FontSize',18)
ylabel('Viral load (log_{10}(copies/mL))')
legend([l1 l2],'Viral load','V(t)')
xlim([0 p.max_time])
xlabel('Days')

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
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)
l2 = plot(t_av,sol_av,'k-','LineWidth',2)
plot([0 p.max_time], [2 2],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'xtick',[0 p.max_time/2 p.max_time],'xticklabel',{'S','S+10','S+20'})
set(gca,'FontSize',18)
ylabel('Viral load (log_{10}(copies/mL))')
legend([l1 l2],'Viral load','V(t)')
xlim([0 p.max_time])
colormap(colmap)
c = colorbar
set(c, 'ticks',linspace(0,1,9),'ticklabels',{'S5','S6','S18','G1','G2','G5','G6','G7',''})
xlabel('Days')

figure
hold on 
yyaxis left
l1 = plot(t_av, sol_S*1e9,'Color',[94,60,153]/255,'LineWidth',2)
l2 = plot(t_av, sol_I*1e9,'Color',[94,60,153]/255,'LineWidth',2)
set(gca,'yscale','log')
ylim([1e0 1e9])
set(gca,'ytick',[1e0 1e5 1e9])
ylabel('Cells/ml')
yyaxis right
l3 = plot(t_av, sol_D*1e9,'Color',[230,97,1]/255,'LineWidth',2)
set(gca,'yscale','log')
ylabel('Cells/ml')
yyaxis right
legend([l1 l2 l3],{'S(t)','I(t)','D(t)'})
set(gca,'yscale','log')
set(gca,'Fontsize',18)
xlabel('Days')
ax = gca;
ax.YAxis(1).Color = [94,60,153]/255;
ax.YAxis(2).Color = [230,97,1]/255;
xlim([0 p.max_time])
set(gca,'xtick',[0 p.max_time/2 p.max_time],'xticklabel',{'S','S+10','S+20'})




%% modifying data to allign time and viral load measurements on specific integer days

time_s5mod =  [1 2 7 8 10 12 13 14 15 16 17 18 19 20 21 22]+2; %adding rounded lags
time_s6mod =  [1 2 4 5 6 7 9 10 11 12 13 14 15 16]+3;
time_s18mod = [1 2 3 4 5 6 7 8 9 10 11 13 14 15 17 19 20]+1;
time_g1mod =  [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19]+1;
time_g2mod =  [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17]+1;
time_g5mod =  [1 2 3 4 5 6 7 8 9 10 11 12 13 14 16]+1;
time_g6mod =  [1 2 3 4 5 6 7 9 10 13 14]+3;
time_g7mod =  [1 2 3 4 7 8 9 10 11 12]+2;

for ii = 1:24

    vecvec = [viral_load_s5(find(time_s5mod == ii));...
        viral_load_s6(find(time_s6mod == ii));...
        viral_load_s18(find(time_s18mod == ii));...
        viral_load_g1(find(time_g1mod == ii));...
        viral_load_g2(find(time_g2mod == ii));...
        viral_load_g5(find(time_g5mod == ii));...
        viral_load_g6(find(time_g6mod == ii));...
        viral_load_g7(find(time_g7mod == ii))];
        
    max_vl(ii) = max(vecvec);
    min_vl(ii) = min(vecvec);
end

figure
hold on
l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[140,81,10]/255,'LineWidth',2.5)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)

plot(1:24,max_vl)
plot(1:24,min_vl)

%% fitting dV dI for max patient and min patient


p.beta = betafit;
p.d_V = d_Vfit;
p.d_I = d_Ifit;
p.V0 = V0fit;

time = 1:24;

datamax = max_vl;

[t_avMAX, sol_avMAX, sol_SMAX, sol_IMAX, sol_DMAX, param_fitMAX] = fit_SIVD_parameters_single(p,datamax,time);

datamin = min_vl;

[t_avMIN, sol_avMIN, sol_SMIN, sol_IMIN, sol_DMIN, param_fitMIN] = fit_SIVD_parameters_single(p,datamin,time);

figure
hold on
l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[140,81,10]/255,'LineWidth',2.5)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)
plot(t_avMAX,sol_avMAX)
plot(t_avMIN,sol_avMIN)
plot(1:24,max_vl)
plot(1:24,min_vl)

SD = abs(sqrt(1)*([param_fitMAX-param_fitMIN])/3.92);

%simulate model for +- SD
p.beta = betafit+SD(1);
p.d_V = d_Vfit+SD(2);
p.d_I = d_Ifit+SD(3);
p.V0 = V0fit+SD(4);

[t_avPLUSSD, sol_avPLUSSD, sol_SPLUSSD, sol_IPLUSSD, sol_DPLUS,solPLUSSD] = simulation_SIVD(p,[0,time])

p.beta = betafit-SD(1);
p.d_V = d_Vfit-SD(2);
p.d_I = d_Ifit-SD(3);
p.V0 = V0fit-SD(4);

[t_avMINUSSD, sol_avMINUSSD, sol_SMINUSSD, sol_IMINSSD, sol_DMINS,solMINSSD] = simulation_SIVD(p,[0,time])

figure
hold on  

LOCS_ABOVE_PLUS = find(t_avPLUSSD>2.9);
LOCS_ABOVE_MINUS = find(t_avMINUSSD>2.9);

inBetween = [[sol_avPLUSSD(1:LOCS_ABOVE_PLUS-1),sol_avMINUSSD(LOCS_ABOVE_PLUS:end)], fliplr([sol_avMINUSSD(1:LOCS_ABOVE_MINUS-1),sol_avPLUSSD(LOCS_ABOVE_PLUS:end)])];
x2 = [t_avPLUSSD, fliplr(t_avMINUSSD)];
fill(x2, inBetween, [1 0.85 0.744],'LineStyle','none');
hold on

l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[0.5 0.5 0.5],'LineWidth',1)
l2 = plot(t_av,sol_av,'Color',[230,97,1]/255,'LineWidth',2)
plot([0 p.max_time], [2 2],'k--','LineWidth',2)
set(gca,'xtick',[0 p.max_time/2,p.max_time],'xticklabel',{'S','S+10','S+20'})
set(gca,'FontSize',18)
ylabel('Viral load (log_{10}(copies/mL))')
legend([l1 l2],'Viral load','V(t)')
xlim([0 p.max_time])
xlabel('Days')

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
inBetween = [[sol_avPLUSSD(1:LOCS_ABOVE_PLUS-1),sol_avMINUSSD(LOCS_ABOVE_PLUS:end)], fliplr([sol_avMINUSSD(1:LOCS_ABOVE_MINUS-1),sol_avPLUSSD(LOCS_ABOVE_PLUS:end)])];
x2 = [t_avPLUSSD, fliplr(t_avMINUSSD)];
fill(x2, inBetween, [0.88 0.88 0.88],'LineStyle','none');
hold on

l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[140,81,10]/255,'LineWidth',2.5)
plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)
l2 = plot(t_av,sol_av,'k-','LineWidth',2)
plot([0 p.max_time], [2 2],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'xtick',[0 p.max_time/2 p.max_time],'xticklabel',{'S','S+10','S+20'})
set(gca,'FontSize',18)
ylabel('Viral load (log_{10}(copies/mL))')
legend([l1 l2],'Viral load','V(t)')
xlim([0 p.max_time])
colormap(colmap)
c = colorbar
set(c, 'ticks',linspace(0,1,9),'ticklabels',{'S5','S6','S18','G1','G2','G5','G6','G7',''})
xlabel('Days')


figure
hold on 
inBetween = [[sol_avPLUSSD(1:LOCS_ABOVE_PLUS-1),sol_avMINUSSD(LOCS_ABOVE_PLUS:end)], fliplr([sol_avMINUSSD(1:LOCS_ABOVE_MINUS-1),sol_avPLUSSD(LOCS_ABOVE_PLUS:end)])];
x2 = [t_avPLUSSD, fliplr(t_avMINUSSD)];
fill(x2, inBetween, [0.88 0.88 0.88],'LineStyle','none');
hold on

l1 = plot(p.lag_s5+time_s5,viral_load_s5,'o','Color',[140,81,10]/255,'LineWidth',2.5)
l2 = plot(p.lag_s6+time_s6,viral_load_s6,'o','Color',[191,129,45]/255,'LineWidth',2.5)
l3 = plot(p.lag_s18+time_s18,viral_load_s18,'o','Color',[223,194,125]/255,'LineWidth',2.5)
l4 = plot(p.lag_g1+time_g1,viral_load_g1,'o','Color',[246,232,195]/255,'LineWidth',2.5)
l5 = plot(p.lag_g2+time_g2,viral_load_g2,'o','Color',[199,234,229]/255,'LineWidth',2.5)
l6 = plot(p.lag_g5+time_g5,viral_load_g5,'o','Color',[128,205,193]/255,'LineWidth',2.5)
l7 = plot(p.lag_g6+time_g6,viral_load_g6,'o','Color',[53,151,143]/255,'LineWidth',2.5)
l8 = plot(p.lag_g7+time_g7,viral_load_g7,'o','Color',[1,102,94]/255,'LineWidth',2.5)
l9 = plot(t_av,sol_av,'k-','LineWidth',2)
plot([0 p.max_time], [2 2],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'xtick',[0 p.max_time/2 p.max_time],'xticklabel',{'S','S+10','S+20'})
set(gca,'FontSize',18)
ylabel('Viral load (log_{10}(copies/mL))')
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9],{'S5','S6','S18','G1','G2','G5','G6','G7','V(t)'})
xlim([0 p.max_time])
%colormap(colmap)
%c = colorbar
%set(c, 'ticks',linspace(0,1,9),'ticklabels',{'S5','S6','S18','G1','G2','G5','G6','G7',''})
xlabel('Days')

