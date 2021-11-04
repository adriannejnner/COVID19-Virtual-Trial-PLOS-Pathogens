function plot_monos_neuts_macs(time,solMMACS,solSMACS,solRMACS,solMMONOS,solSMONOS,solRMONOS,solMNEUTS,solSNEUTS,solRNEUTS,p)

% knockdown plots -----------------------------------------------------------

figure
hold on 
plot(time, solMNEUTS(1,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solRMONOS(1,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solRMACS(1,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solRNEUTS(1,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Viral load')
ylabel('log_{10}(copies/ml)')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7A_KD.fig');   
saveas(gcf,'Fig_7A_KD.png');


figure
hold on 
plot(time, (solMNEUTS(2,:)+solMMONOS(4,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solRMONOS(2,:)+solRMONOS(4,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solRMACS(2,:)+solRMACS(4,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solRNEUTS(2,:)+solRNEUTS(4,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','log')
title('Uninfected cells')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7B_KD.fig');   
saveas(gcf,'Fig_7B_KD.png');


figure
hold on 
plot(time, (solMNEUTS(7,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solRMONOS(7,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solRMACS(7,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solRNEUTS(7,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('inflammatory macrophages')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7C_KD.fig');   
saveas(gcf,'Fig_7C_KD.png');

figure
hold on 
plot(time, (solMNEUTS(9,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solRMONOS(9,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solRMACS(9,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solRNEUTS(9,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Neutrophils')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7D_KD.fig');   
saveas(gcf,'Fig_7D_KD.png');



figure
hold on 
plot(time, (solMNEUTS(10,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solRMONOS(10,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solRMACS(10,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solRNEUTS(10,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cells')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_1_KD.fig');   
saveas(gcf,'Fig_7E_1_KD.png');

figure
hold on 
plot(time, solMNEUTS(10,:)./solMMONOS(3,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solRMONOS(10,:)./solRMONOS(3,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solRMACS(10,:)./solRMACS(3,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solRNEUTS(10,:)./solRNEUTS(3,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cells:uninfected cells')
ylabel('ratio')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_2_KD.fig');   
saveas(gcf,'Fig_7E_2_KD.png');


figure
hold on 
plot(time, solMNEUTS(10,:)./solMMONOS(11,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solRMONOS(10,:)./solMMONOS(11,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solRMACS(10,:)./solMMONOS(11,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solRNEUTS(10,:)./solMMONOS(11,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cell:IL-6')
ylabel('ratio')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_3_KD.fig');   
saveas(gcf,'Fig_7E_3_KD.png');


figure
hold on 
plot(time, solMNEUTS(11,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solRMONOS(11,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solRMACS(11,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solRNEUTS(11,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Unbound IL-6')
ylabel('pg/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7F_KD.fig');   
saveas(gcf,'Fig_7F_KD.png');

% knockdown plots -----------------------------------------------------------


figure
hold on 
plot(time, solMNEUTS(1,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solSMONOS(1,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solSMACS(1,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solSNEUTS(1,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Viral load')
ylabel('log_{10}(copies/ml)')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7A_KO.fig');   
saveas(gcf,'Fig_7A_KO.png');


figure
hold on 
plot(time, (solMNEUTS(2,:)+solMMONOS(4,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solSMONOS(2,:)+solRMONOS(4,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solSMACS(2,:)+solRMACS(4,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solSNEUTS(2,:)+solRNEUTS(4,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','log')
title('Uninfected cells')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7B_KO.fig');   
saveas(gcf,'Fig_7B_KO.png');


figure
hold on 
plot(time, (solMNEUTS(7,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solSMONOS(7,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solSMACS(7,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solSNEUTS(7,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('cells/ml')
ylabel('inflammatory macrophages')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7C_KO.fig');   
saveas(gcf,'Fig_7C_KO.png');


figure
hold on 
plot(time, (solMNEUTS(9,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solSMONOS(9,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solSMACS(9,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solSNEUTS(9,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Neutrophils')
ylabel('cells/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7D_KO.fig');   
saveas(gcf,'Fig_7D_KO.png');



figure
hold on 
plot(time, (solMNEUTS(10,:))*1e9,'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, (solSMONOS(10,:))*1e9,'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, (solSMACS(10,:))*1e9,':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, (solSNEUTS(10,:))*1e9,'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cells')
ylabel('cells/ml')
set(gca,'Fontsize',24)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_1_KO.fig');   
saveas(gcf,'Fig_7E_1_KO.png');

figure
hold on 
plot(time, solMNEUTS(10,:)./solMMONOS(3,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solSMONOS(10,:)./solMMONOS(3,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solSMACS(10,:)./solMMONOS(3,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solSNEUTS(10,:)./solMMONOS(3,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cells:infected cells')
ylabel('ratio')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_2_KO.fig');   
saveas(gcf,'Fig_7E_2_KO.png');


figure
hold on 
plot(time, solMNEUTS(10,:)*1e9./solMMONOS(11,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solSMONOS(10,:)*1e9./solMMONOS(11,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solSMACS(10,:)*1e9./solMMONOS(11,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solSNEUTS(10,:)*1e9./solMMONOS(11,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('CD8^+ T cells/IL-6')
ylabel('ratio')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7E_3_KO.fig');   
saveas(gcf,'Fig_7E_3_KO.png');


figure
hold on 
plot(time, solMNEUTS(11,:),'Color',[108 191 208]/255,'LineWidth',3) % original parameters simulation
plot(time, solSMONOS(11,:),'-.','Color',[133 63 63]/255,'LineWidth',3);  % monocytes
plot(time, solSMACS(11,:),':','Color',[242,180,180]/255,'LineWidth',3); % macrophages 
plot(time, solSNEUTS(11,:),'--','Color',[212 102 131]/255,'LineWidth',3);  % neutrophils
set(gca,'yscale','linear')
title('Unbound IL-6')
ylabel('pg/ml')
set(gca,'Fontsize',22)
xlabel('Time (days)')
saveas(gcf,'Fig_7F_KO.fig');   
saveas(gcf,'Fig_7F_KO.png');

end