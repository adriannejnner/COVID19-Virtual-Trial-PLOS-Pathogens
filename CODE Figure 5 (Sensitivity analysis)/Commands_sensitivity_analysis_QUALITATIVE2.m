
tol = 5;
tspan = [0 30];
[timeORIG,solORIG] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model

factor = 0.2;% 0.1 - 10%

col = [165,0,38;...
215,48,39;...
244,109,67;...
253,174,97;...
254,224,144;...
255,255,191;...
224,243,248;...
171,217,233;...
116,173,209;...
69,117,180;...
49,54,149]/255;%jet(11);

p.del_V_N = p.del_V_N*(1-factor);
[time1,sol1] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.del_V_N = p.del_V_N/(1-factor);

p.p_F_I = p.p_F_I*(1-factor);
[time2,sol2] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_F_I = p.p_F_I/(1-factor);

p.beta = p.beta*(1-factor);
[time3,sol3] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.beta = p.beta/(1-factor);

p.phat = p.phat*(1-factor);
[time4,sol4] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.phat = p.phat/(1-factor);

figure
hold on 
plot(time1,sol1(1,:),'Color',col(11,:),'LineWidth',3)
plot(time2,sol2(1,:),'Color',col(8,:),'LineWidth',3)
plot(time3,sol3(1,:),'Color',col(5,:),'LineWidth',3)
plot(time4,sol4(1,:),'Color',col(1,:),'LineWidth',3)
legend('\delta_{V,N}','p_{F,M\Phi}','\beta','p')
xlabel('Time (days)')
ylabel('log_{10}(copies/ml)')
title('Viral load')
set(gca,'FontSize',22)
axes('position',[0.4 .575 .25 .25])
box on
hold on
plot(time1,sol1(1,:),'Color',col(11,:),'LineWidth',3)
plot(time2,sol2(1,:),'Color',col(8,:),'LineWidth',3)
plot(time3,sol3(1,:),'Color',col(5,:),'LineWidth',3)
plot(time4,sol4(1,:),'Color',col(1,:),'LineWidth',3)
xlim([1,3])
axis tight

%%

p.k_lin_F = p.k_lin_F*(1-factor);
[time1,sol1] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.k_lin_F = p.k_lin_F/(1-factor);

p.p_F_I = p.p_F_I*(1-factor);
[time2,sol2] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_F_I = p.p_F_I/(1-factor);

p.p_F_MPhi = p.p_F_MPhi*(1-factor);
[time3,sol3] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_F_MPhi = p.p_F_MPhi/(1-factor);

p.phat = p.phat*(1-factor);
[time4,sol4] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.phat = p.phat/(1-factor);

figure
hold on 
plot(time1,(sol1(2,:)+sol1(4,:))*1e9,'Color',col(11,:),'LineWidth',3)
plot(time4,(sol4(2,:)+sol4(4,:))*1e9,'Color',col(7,:),'LineWidth',3)
plot(time3,(sol3(2,:)+sol3(4,:))*1e9,'Color',col(5,:),'LineWidth',3)
plot(time2,(sol2(2,:)+sol2(4,:))*1e9,'Color',col(3,:),'LineWidth',3)
legend('k_{lin}_F','p','p_{F,M\Phi}','p_{F,I}')
xlabel('Time (days)')
ylabel('cells/ml')
title('Uninfected cells')
set(gca,'FontSize',22)
set(gca,'yscale','log')
axes('position',[0.4 .28 .25 .25])
box on
hold on
plot(time1,(sol1(2,:)+sol1(4,:))*1e9,'Color',col(11,:),'LineWidth',3)
plot(time4,(sol4(2,:)+sol4(4,:))*1e9,'Color',col(7,:),'LineWidth',3)
plot(time3,(sol3(2,:)+sol3(4,:))*1e9,'Color',col(5,:),'LineWidth',3)
plot(time2,(sol2(2,:)+sol2(4,:))*1e9,'Color',col(3,:),'LineWidth',3)
xlim([1.5,3.5])
axis tight

%%

p.k_int_L = p.k_int_L*(1-factor);
[time1,sol1] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.k_int_L = p.k_int_L/(1-factor);

p.k_B_L = p.k_B_L*(1-factor);
[time2,sol2] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.k_B_L = p.k_B_L/(1-factor);

p.p_L_MPhi = p.p_L_MPhi*(1-factor);
[time3,sol3] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_L_MPhi = p.p_L_MPhi/(1-factor);

p.del_V_N = p.del_V_N*(1-factor);
[time4,sol4] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.del_V_N = p.del_V_N/(1-factor);

figure
hold on 
plot(time2,sol2(11,:),'Color',col(11,:),'LineWidth',3)
plot(time4,sol4(11,:),'Color',col(7,:),'LineWidth',3)
plot(time3,sol3(11,:),'Color',col(4,:),'LineWidth',3)
plot(time1,sol1(11,:),'Color',col(2,:),'LineWidth',3)
legend('k_B_L','\delta_{V,N}','p_{L,M\Phi}','k_{lin_L}')
xlabel('Time (days)')
ylabel('pg/ml')
title('Unbound IL-6')
set(gca,'FontSize',22)



%%


p.k_lin_F = p.k_lin_F*(1-factor);
[time1,sol1] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.k_lin_F = p.k_lin_F/(1-factor);

p.p_M_I = p.p_M_I*(1-factor);
[time2,sol2] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_M_I = p.p_M_I/(1-factor);

p.del_I_T = p.del_I_T*(1-factor);
[time3,sol3] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.del_I_T = p.del_I_T/(1-factor);

p.p_F_I = p.p_F_I*(1-factor);
[time4,sol4] = COVID_IMMUNE_MODEL_COPY(p,tspan);% solve model
p.p_F_I = p.p_F_I/(1-factor);

figure
hold on 
plot(time1,sol1(17,:),'Color',col(11,:),'LineWidth',3)
plot(time3,sol3(17,:),'Color',col(9,:),'LineWidth',3)
plot(time2,sol2(17,:),'Color',col(4,:),'LineWidth',3)
plot(time4,sol4(17,:),'Color',col(2,:),'LineWidth',3)
legend('k_{lin_F}','\delta_{I,T}','p_{M,I}','p_{F,I}')
xlabel('Time (days)')
ylabel('pg/ml')
title('Unbound IFN')
set(gca,'FontSize',22)
