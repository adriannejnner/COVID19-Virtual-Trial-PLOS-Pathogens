
load('patients_sim_annealing_severeNEW2.mat') % if no good, try 4 or no number
patient = patients_new;
mu_vec =    [p.beta,    p.p_MPhi_I_L,    p.p_L_MPhi,    p.p_F_I,    p.eta_F_I,    p.eps_L_T,    p.p_M_I,    p.eta_F_MPhi,    p.tau_T     p.eps_F_I,    p.p_F_M];

figure
hold on 
histogram(patient(:,2),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,2)) mean(patient(:,2))],[0 60],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(2) mu_vec(2)],[0 60],'k:','LineWidth',3)
xlabel('p_{M\Phi_I,L}')
set(gca,'FontSize',18)
title({'Rate of monocyte to macrophage','differentiation by IL-6'})
saveas(gcf,'Fig_S10A.fig');
saveas(gcf,'Fig_S10A.png');

figure
hold on 
histogram(patient(:,3),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,3)) mean(patient(:,3))],[0 50],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(3) mu_vec(3)],[0 50],'k:','LineWidth',3)
xlabel('p_{L,M\Phi}')
set(gca,'FontSize',18)
title({'Rate of IL-6 production by','activated macrophages'})
saveas(gcf,'Fig_S10B.fig');
saveas(gcf,'Fig_S10B.png');

figure
hold on 
histogram(patient(:,4),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,4)) mean(patient(:,4))],[0 50],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(4) mu_vec(4)],[0 50],'k:','LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{F,I}')
title({'Production rate of IFN','by infected cells'})
saveas(gcf,'Fig_S10C.fig');
saveas(gcf,'Fig_S10C.png');

figure
hold on 
histogram(patient(:,7),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,7)) mean(patient(:,7))],[0 60],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(7) mu_vec(7)],[0 60],'k:','LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{M,I}')
title({'Monocyte recruitment rate','by infected cells'})
saveas(gcf,'Fig_S10D.fig');
saveas(gcf,'Fig_S10D.png');

figure
hold on 
histogram(patient(:,8),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,8)) mean(patient(:,8))],[0 100],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(8) mu_vec(8)],[0 100],'k:','LineWidth',3)
set(gca,'FontSize',18)
xlabel('\eta_{F,M\Phi}')
title({'IFN production by inflammatory','macrophages half-effect'})
saveas(gcf,'Fig_S10E.fig');
saveas(gcf,'Fig_S10E.png');

figure
hold on 
histogram(patient(:,10),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,10)) mean(patient(:,10))],[0 60],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(10) mu_vec(10)],[0 60],'k:','LineWidth',3)
set(gca,'FontSize',18)
xlabel('\epsilon_{F,I}')
title({'IFN inhibition of viral production','half-effect'})
saveas(gcf,'Fig_S10F.fig');
saveas(gcf,'Fig_S10F.png');

figure
hold on 
histogram(patient(:,11),'FaceColor',[226,217,243]/255,'EdgeColor',[169,120,186]/255)
plot([mean(patient(:,11)) mean(patient(:,11))],[0 60],'--','Color',[0.49,0.18,0.56],'LineWidth',3)
plot([mu_vec(11) mu_vec(11)],[0 60],'k:','LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{M,F}')
title({'Monocyte production rate','of IFN'})
saveas(gcf,'Fig_S10G.fig');
saveas(gcf,'Fig_S10G.png');
