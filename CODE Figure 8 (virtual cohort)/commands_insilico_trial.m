
tspan=[0 21];

%p = load_parameters_simulated_annealing2();

%p = Homeostasis_calculations(p);

patients_full_cohort = [];
load('patients_sim_annealing_severe12.mat')%2 , %6 is the right one from the original submission
patients_full_cohort(1:3,:) = patient([306,391,400],:);
load('patients_sim_annealing_severe13.mat')%2 , %6 is the right one from the original submission
patients_full_cohort = [patients_full_cohort;patient];%([12    86   184   187   200   343   368]-4,:)];
load('patients_sim_annealing_severe14.mat')%2 , %6 is the right one from the original submission
patients_full_cohort = [patients_full_cohort;patient];%([12    86   184   187   200   343   368]-4,:)];
load('patients_sim_annealing_severe15.mat')%2 , %6 is the right one from the original submission
patients_full_cohort = [patients_full_cohort;patient];%([12    86   184   187   200   343   368]-4,:)];

%load('patients_sim_annealing_MC.mat')
%patients_full_cohort = [patients_full_cohort;patient];

time_grid = linspace(tspan(1),tspan(2),1000);

p.p_F_IORIG = p.p_F_I;
p.eta_F_IORIG = p.eta_F_I;
p.p_M_IORIG = p.p_M_I;
p.eta_F_MPhiORIG = p.eta_F_MPhi;

p.I0 = 0;
%p.Vtrans_0 = 5.21;
%p.V0 = p.Vtrans_0;

%p.beta = 0.8;
%p.d_I = 5.19;
%p.p = 72.94;
%p.d_V = 0.07;
%p.ktr =30.32;

for jj = 1:length(patients_full_cohort);%   *15  *26 *10 *10
    
   % p.beta = patients_full_cohort(jj,1);
    p.p_MPhi_I_L = patients_full_cohort(jj,2);
    p.p_L_MPhi = patients_full_cohort(jj,3);
    p.p_F_I = patients_full_cohort(jj,4);
    %p.eta_F_I = patients_full_cohort(jj,5);
    %p.eps_L_T = patients_full_cohort(jj,6);
    p.p_M_I = patients_full_cohort(jj,7);
    p.eta_F_MPhi = patients_full_cohort(jj,8);
   % p.tau_T = patients_full_cohort(jj,9);
    p.eps_F_I =  patients_full_cohort(jj,10);
    p.p_F_M = patients_full_cohort(jj,11);

    %p = Homeostasis_calculations(p);
    
    [time,sol,solstruc] = COVID_IMMUNE_MODELINSILICO(p,tspan);
   
    sol_virus(jj,:) = deval(solstruc,time_grid,1);
    sol_tissue(jj,:) = deval(solstruc,time_grid,2)+deval(solstruc,time_grid,4);
    sol_infected(jj,:) = deval(solstruc,time_grid,3);
    sol_dead(jj,:) = deval(solstruc,time_grid,5);
    sol_macs_res(jj,:) = deval(solstruc,time_grid,6);
    sol_macs_inflam(jj,:) = deval(solstruc,time_grid,7);
    sol_monocytes(jj,:) = deval(solstruc,time_grid,8);
    sol_neutrophils(jj,:) = deval(solstruc,time_grid,9);
    sol_Tcells(jj,:) = deval(solstruc,time_grid,10);
    
    sol_IL6(jj,:) = deval(solstruc,time_grid,11); 
    sol_IFN(jj,:) = deval(solstruc,time_grid,17); 
    sol_GCSF(jj,:) = deval(solstruc,time_grid,15);
    
    days14 = find(time>14,1);
    
    if isreal(sol)==1
        minimum_tissue(jj) = min(sol(2,:)+sol(4,:));%+sol(3,:)
        
        max_IL6(jj) = max(sol(11,:));
        max_bound_IL6(jj) = max(sol(12,:));
        max_dead_tissue(jj) = max(sol(5,:));
        max_inflam_macs(jj) = max(sol(7,:));

        tgrid_vec = [];
        for i = 1:length(sol(17,:))-1
           tgrid_vec(i) = time(i+1)-time(i); 
        end
        IFN_exposure(jj) = sum(sol(17,1:end-1).*tgrid_vec);

        peak_loc = find(sol(17,:)==max(sol(17,:)));
        IFN_peak(jj) = time(peak_loc); %time of IFN peak
        max_T_cells(jj) = max(sol(10,1:days14));    
        max_infected_cells(jj) = max(sol(3,:));
        max_inflam_macs(jj) = max(sol(7,:));
        max_neutrophils(jj) = max(sol(9,:));
        max_monocytes(jj) = max(sol(8,:));
    else
        minimum_tissue(jj) = NaN;%+sol(3,:)
        max_IL6(jj) = NaN;
        max_bound_IL6(jj) = NaN;
        max_dead_tissue(jj) = NaN;
        IFN_exposure(jj) = NaN;
        IFN_peak(jj) = NaN; %time of IFN peak
        max_T_cells(jj) = NaN;    
        max_infected_cells(jj) = NaN;
        max_inflam_macs(jj) = NaN;
        max_neutrophils(jj) = NaN;
        max_monocytes(jj) = NaN;
    end
    jj
end

save('Patient responsesALL.mat','time_grid','minimum_tissue','max_IL6','IFN_exposure','IFN_peak','max_T_cells','max_inflam_macs','max_monocytes','sol_virus','sol_tissue','sol_infected','sol_dead','sol_macs_res','sol_macs_inflam','sol_neutrophils','sol_monocytes','sol_Tcells','sol_IL6','sol_IFN','sol_GCSF')
find(IFN_peak>11)

locs_IFNpeak = find(IFN_peak>3.8 & IFN_peak <10.5);
patients_withoutIFN = patients_full_cohort;
patients_withoutIFN(locs_IFNpeak(1:70),:) = [];

remaining_numberpatients = 200-length(locs_IFNpeak(1:70));

%rng(0,'twister');
r = randi([1 length(patients_withoutIFN)],1,remaining_numberpatients);

patients_new = [patients_full_cohort(locs_IFNpeak(1:70),:);patients_withoutIFN(r,:)];

 save('patients_sim_annealing_severeNEW.mat','patients_new')
 
 minimum_tissueNEW = [minimum_tissue(locs_IFNpeak(1:70)),minimum_tissue(r)];
 max_IL6NEW = [ max_IL6(locs_IFNpeak(1:70)), max_IL6(r)];
 max_neutrophilsNEW = [max_neutrophils(locs_IFNpeak(1:70)),max_neutrophils(r)];
 max_inflam_macsNEW = [max_inflam_macs(locs_IFNpeak(1:70)),max_inflam_macs(r)];
 IFN_peakNEW = [IFN_peak(locs_IFNpeak(1:70)),IFN_peak(r)];
 max_T_cellsNEW = [max_T_cells(locs_IFNpeak(1:70)),max_T_cells(r)];
 
 PsiNEW = (max_IL6NEW)/(mean(max_IL6NEW))+(max_neutrophilsNEW)/mean(max_neutrophilsNEW)+(p.Smax-minimum_tissueNEW)/mean(p.Smax-minimum_tissueNEW);
[sorted,where] = sort(real(PsiNEW));
%% Ordering patients by surragate inflammation marker

%load('Patient responsesALL.mat')
tspan=[0 21];
time_grid = linspace(tspan(1),tspan(2),1000);
p.Smax = 0.16;
%[max_IL6 max_neutrophils minimum_tissue, max_inflam_macs, IFN_exposure, IFN_peak, max_T_cells]=reevaluate_patients(time_grid,sol_virus,sol_tissue,sol_infected,sol_dead,sol_macs_res,sol_macs_inflam,sol_neutrophils, sol_monocytes,sol_Tcells, sol_IL6, sol_IFN, sol_GCSF,time_grid);

%surragate marker for severity 
Psi = (max_IL6)/(mean(max_IL6))+(max_neutrophils)/mean(max_neutrophils)+(p.Smax-minimum_tissue)/mean(p.Smax-minimum_tissue);
[sorted,where] = sort(real(Psi));

figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(real(max_IL6(where)),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('max IL-6 (pg/ml)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,real(max_IL6(where)));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8A.fig');
saveas(gcf,'Fig_8A.png');


figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(real(max_neutrophils(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)%max_bound_IL6(where)),'o:')
ylabel('max neutrophils (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,real(max_neutrophils(where)));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8B.fig');
saveas(gcf,'Fig_8B.png');

figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(real(max_inflam_macs(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('max inflam. macs (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,real(max_inflam_macs(where)));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8C.fig');
saveas(gcf,'Fig_8C.png');

figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(IFN_exposure(where),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('IFN exposure (pg/ml)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,IFN_exposure(where));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8D.fig');
saveas(gcf,'Fig_8D.png');

figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(real(IFN_peak(where)),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('IFN peak (days)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,IFN_peak(where));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8E.fig');
saveas(gcf,'Fig_8E.png');

figure
yyaxis left
hold on
plot(sorted,'o','Color',[127,205,187]/255,'LineWidth',1)
ylabel('Inflammation marker, \Psi')
yyaxis right
plot(real(max_T_cells(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('Max T cells (cells/ml)')
xlabel('Patients (ordered by \Psi)')
set(gca,'FontSize',18)
ax = gca;
ax.YAxis(1).Color = [0.33 0.66 0.58];
ax.YAxis(2).Color = [44,127,184]/255;
[R,P] = corrcoef(sorted,real(max_T_cells(where)));
yyaxis left
text(2,3.85,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
yyaxis right
ylh = get(gca,'ylabel');
set(ylh,'Rotation',270);
saveas(gcf,'Fig_8F.fig');
saveas(gcf,'Fig_8F.png');

%-------------------------------------------------

figure
hold on
plot(sorted,real(max_IL6(where)),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('max IL-6 (pg/ml)')
xlabel('Patient Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_IL6(where)));
text(2.05,93,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8A_2.fig');
saveas(gcf,'Fig_8A_2.png');

figure
hold on
plot(sorted,real(max_neutrophils(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)%max_bound_IL6(where)),'o:')
ylabel('max neutrophils (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_neutrophils(where)));
text(2.05,10.55e6,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8B_2.fig');
saveas(gcf,'Fig_8B_2.png');

figure
hold on
plot(sorted,real(max_inflam_macs(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('max inflam. macs (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_inflam_macs(where)));
text(2.05,2.75*1e5,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8C_2.fig');
saveas(gcf,'Fig_8C_2.png');

figure
hold on
plot(sorted,IFN_exposure(where),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('IFN exposure (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,IFN_exposure(where));
text(2.05,5.45,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8D_2.fig');
saveas(gcf,'Fig_8D_2.png');

figure
hold on
plot(sorted,real(IFN_peak(where)),'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('IFN peak (days)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,IFN_peak(where));
text(2.05,13.5,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8E_2.fig');
saveas(gcf,'Fig_8E_2.png');

figure
hold on
plot(sorted,real(max_T_cells(where))*1e9,'o','Color',[44,127,184]/255,'LineWidth',1)
ylabel('Max T cells (cells/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_T_cells(where)));
text(2.05,3.7*1e6,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8F_2.fig');
saveas(gcf,'Fig_8F_2.png');

%-------------------------------------------------

colvec = [237,248,177;...
199,233,180;...
127,205,187;...
65,182,196;...
29,145,192;...
34,94,168;...
37,52,148;...
8,29,88]/255;

colgrid = linspace(min(sorted),max(sorted),9);

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),real(max_IL6(where(ii))),'o','Color',col,'LineWidth',2)
end
ylabel('max IL-6 (pg/ml)')
xlabel('Patient Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_IL6(where)));
text(2.05,93,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8A_3.fig');
saveas(gcf,'Fig_8A_3.png');

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),real(max_neutrophils(where(ii)))*1e9,'o','Color',col,'LineWidth',2)%max_bound_IL6(where)),'o:')
end
    ylabel('max neutrophils (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_neutrophils(where)));
text(2.05,10.55e6,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8B_3.fig');
saveas(gcf,'Fig_8B_3.png');

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),real(max_inflam_macs(where(ii)))*1e9,'o','Color',col,'LineWidth',2)
end
ylabel('max inflam. macs (cells/ml)')%max IL-6 bound concentration (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_inflam_macs(where)));
text(2.05,2.75*1e5,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8C_3.fig');
saveas(gcf,'Fig_8C_3.png');

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),IFN_exposure(where(ii)),'o','Color',col,'LineWidth',2)
end
ylabel('IFN exposure (pg/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,IFN_exposure(where));
text(2.05,5.45,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8D_3.fig');
saveas(gcf,'Fig_8D_3.png');

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),real(IFN_peak(where(ii))),'o','Color',col,'LineWidth',2)
end
ylabel('IFN peak (days)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,IFN_peak(where));
text(2.05,13.5,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8E_3.fig');
saveas(gcf,'Fig_8E_3.png');

figure
hold on
for ii = 1:603;%size(patients_full_cohort,1)
    loc = find(colgrid>sorted(ii),1);
    if(isempty(loc)==1)
        col = colvec(end,:);
    else
        col = colvec(loc-1,:);
    end
    plot(sorted(ii),real(max_T_cells(where(ii)))*1e9,'o','Color',col,'LineWidth',2)
end
ylabel('Max T cells (cells/ml)')
xlabel('Patients Inflammation marker \Psi')
set(gca,'FontSize',18)
[R,P] = corrcoef(sorted,real(max_T_cells(where)));
text(2.05,3.7*1e6,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
saveas(gcf,'Fig_8F_3.fig');
saveas(gcf,'Fig_8F_3.png');

figure
colormap(colvec)
c = colorbar
set(gca,'FontSize',21)
set(c,'ticks',linspace(0,1,9),'ticklabels',{'2.2','2.4','2.6','2.8','3','3.3','3.5','3.7','3.9'})


%---------------------------------
patient = patients_full_cohort;
colmap = [50,136,189;...
102,194,165;...
171,221,164;...
230,245,152;...
255,255,191;...
254,224,139;...
253,174,97;...
244,109,67;...
213,62,79]/255;

IFN_peak_vec = linspace(min(IFN_peak),max(IFN_peak),10);
figure
hold on 
for ij = 1:603%length(patient)
   
   colval = find(IFN_peak_vec>IFN_peak(where(ij)),1)
    if isempty(colval) ==1
       colval = 9; 
    end
    plot(ij,patient(where(ij),2),'.','Color',colmap(colval-1,:),'MarkerSize',30) 
end
xlabel('Patients ordered by inflammation \Psi')
ylabel('p_{M\Phi_I,L}')
set(gca,'FontSize',18)
[R,P] = corrcoef(IFN_peak(where),patient(where,2));
text(2,14,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
colormap(colmap)
c = colorbar;
set(c,'ticks',linspace(0,1,5),'ticklabels',{'1', '2.4', '3.8', '5.2', '6.6','8','9.5','10.9','12.3','14'})
saveas(gcf,'Fig_8I.fig');
saveas(gcf,'Fig_8I.png');

figure
hold on 
for ij = 1:603%length(patient)
   
   colval = find(IFN_peak_vec>IFN_peak(where(ij)),1)
    if isempty(colval) ==1
       colval = 9; 
    end
    plot(ij,patient(where(ij),4),'.','Color',colmap(colval-1,:),'MarkerSize',30) 
end
xlabel('Patients ordered by inflammation \Psi')
ylabel('p_{F,I}')
set(gca,'FontSize',18)
[R,P] = corrcoef(IFN_peak(where),patient(where,2));
text(2,7.5,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
colormap(colmap)
c = colorbar;
set(c,'ticks',linspace(0,1,10),'ticklabels',{'1', '2.4', '3.8', '5.2', '6.6','8','9.5','10.9','12.3','14'})
saveas(gcf,'Fig_8J.fig');
saveas(gcf,'Fig_8J.png');

max_IFN = max(sol_IFN');
colmap = [77,0,75;...
129,15,124;...
136,65,157;...
140,107,177;...
140,150,198;...
158,188,218;...
191,211,230;...
224,236,244;...
247,252,253]/255;

Tcell_grid = linspace(min(max_T_cells),max(max_T_cells),10);
figure
hold on 
for i = 1:603%length(patient)
    
    colindex = find(Tcell_grid>max_T_cells(where(i)),1)-1;
    if isempty(colindex)==1
    colindex = 9;
    end
    plot(max_IFN(where(i)),max_IL6(where(i)),'.','Color',colmap(colindex,:),'LineWidth',14,'MarkerSize',40)
end

grid on 
colormap(colmap)
c = colorbar;
set(c,'ticks',linspace(0,1,10),'ticklabels',{'0.4','0.7','1','1.3','1.7','2','2.3','2.7','3','3.3'})
set(gca,'FontSize',18)
xlabel('Max IFN (pg/ml)')
ylabel('Max IL-6 (pg/ml)')
saveas(gcf,'Fig_8G.fig');
saveas(gcf,'Fig_8G.png');

%% Violin plots for cohort

load('Patient responses7.mat')
load('Lucas_data.mat')

p.Smax = 0.16;
[max_IL6 max_neutrophils minimum_tissue, max_inflam_macs, IFN_exposure, IFN_peak, max_T_cells, max_monocytes]=reevaluate_patients(time_grid,sol_virus,sol_tissue,sol_infected,sol_dead,sol_macs_res,sol_macs_inflam,sol_neutrophils, sol_monocytes,sol_Tcells, sol_IL6, sol_IFN, sol_GCSF,time_grid);

%surragate marker for severity 
Psi = (max_IL6)/(mean(max_IL6))+(max_neutrophils)/mean(max_neutrophils)+(p.Smax-minimum_tissue)/mean(p.Smax-minimum_tissue);
[sorted,where] = sort(real(Psi));

%violin plots of 50 most severe and 50 most mild (using Psi) AND plotting
%the maximum values for monocytes, neutrophils and IL6
violin_mat_monocyte = [[Monocytes_Lucas_moderate;NaN(2,1)]./mean(Monocytes_Lucas_HCW), max_monocytes(where(1:50))'/p.M0,[Monocytes_Lucas_Severe;NaN(50-29,1)]./mean(Monocytes_Lucas_HCW),max_monocytes(where(end-49:end))'/p.M0];
figure
hold on 
violinplot(violin_mat_monocyte)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Monocytes')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')

violin_mat_neutrophils = [Neutrophils_Lucas_moderate./mean(Neutrophils_Lucas_HCW), [max_neutrophils(where(1:50))'/p.N0;NaN(57-50,1)], [Neutrophils_Lucas_Severe;NaN(57-32,1)]./mean(Neutrophils_Lucas_HCW), [max_neutrophils(where(end-49:end))'/p.N0;NaN(57-50,1)]];
figure
hold on 
violinplot(violin_mat_neutrophils)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Neutrophils')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')

violin_mat_IL6 = [10.^IL6_Lucas_moderate_individual'./mean(10.^IL6_Lucas_HCW_individual),[max_IL6(where(1:50))'/p.L_U_0;NaN(108-50,1)],[10.^IL6_Lucas_Severe_individual'./mean(10.^IL6_Lucas_HCW_individual);NaN(108-41,1)],[max_IL6(where(end-49:end))'/p.L_U_0;NaN(108-50,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('IL-6')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','linear')

%violin plots of 50 most severe and 50 most mild (using Psi) AND plotting
%the maximum values for monocytes, neutrophils and IL6
violin_mat_monocyte = [[Monocytes_Lucas_moderate;NaN(2,1)], max_monocytes(where(1:50))'*1e9,[Monocytes_Lucas_Severe;NaN(50-29,1)],max_monocytes(where(end-49:end))'*1e9];
figure
hold on 
violinplot(violin_mat_monocyte)
ylabel('Cells/ml')
set(gca,'FontSize',18)
title('Monocytes')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')

violin_mat_neutrophils = [Neutrophils_Lucas_moderate, [max_neutrophils(where(1:50))';NaN(57-50,1)]*1e9, [Neutrophils_Lucas_Severe;NaN(57-32,1)], [max_neutrophils(where(end-49:end))';NaN(57-50,1)]*1e9];
figure
hold on 
violinplot(violin_mat_neutrophils)
ylabel('Cells/ml')
set(gca,'FontSize',18)
title('Neutrophils')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


violin_mat_IL6 = [10.^IL6_Lucas_moderate_individual',[max_IL6(where(1:50))';NaN(108-50,1)],[10.^IL6_Lucas_Severe_individual';NaN(108-41,1)],[max_IL6(where(end-49:end))';NaN(108-50,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
ylabel('pg/ml')
set(gca,'FontSize',18)
title('IL-6')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')

day_5 = find(time_grid>=5,1);
monocyte_sample_day5 = sol_monocytes(:,day_5);
%the maximum values for monocytes, neutrophils and IL6
violin_mat_monocyte = [[Monocytes_Lucas_moderate;NaN(2,1)]./mean(Monocytes_Lucas_HCW), monocyte_sample_day5(where(1:50))/p.M0,[Monocytes_Lucas_Severe;NaN(50-29,1)]./mean(Monocytes_Lucas_HCW),monocyte_sample_day5(where(end-49:end))/p.M0];
figure
hold on 
violinplot(violin_mat_monocyte)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Monocytes')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


neutrophil_sample_day5 = sol_neutrophils(:,day_5);
violin_mat_neutrophils = [Neutrophils_Lucas_moderate./mean(Neutrophils_Lucas_HCW), [neutrophil_sample_day5(where(1:50))/p.N0;NaN(57-50,1)], [Neutrophils_Lucas_Severe;NaN(57-32,1)]./mean(Neutrophils_Lucas_HCW), [neutrophil_sample_day5(where(end-49:end))/p.N0;NaN(57-50,1)]];
figure
hold on 
violinplot(violin_mat_neutrophils)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Neutrophils')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


IL6_sample_day5 = sol_IL6(:,day_5);
violin_mat_IL6 = [10.^IL6_Lucas_moderate_individual'./mean(10.^IL6_Lucas_HCW_individual),[IL6_sample_day5(where(1:50))/p.L_U_0;NaN(108-50,1)],[10.^IL6_Lucas_Severe_individual'./mean(10.^IL6_Lucas_HCW_individual);NaN(108-41,1)],[IL6_sample_day5(where(end-49:end))/p.L_U_0;NaN(108-50,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('IL-6')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','linear')

monocyte_sample_dayrand =[];
neutrophil_sample_dayrand = [];
IL6_sample_dayrand = [];

day_5 = find(time_grid>=3,1)
day_15 = find(time_grid>=10,1);
times = randi([day_5 day_15],1,200);
for ii = 1:200
    monocyte_sample_dayrand(ii) = sol_monocytes(ii,times(ii));
    neutrophil_sample_dayrand(ii) = sol_neutrophils(ii,times(ii));
    IL6_sample_dayrand(ii) = sol_IL6(ii,times(ii));
end

%the maximum values for monocytes, neutrophils and IL6
violin_mat_monocyte = [[Monocytes_Lucas_moderate;NaN(2,1)]./mean(Monocytes_Lucas_HCW), monocyte_sample_dayrand(where(1:50))'/p.M0,[Monocytes_Lucas_Severe;NaN(50-29,1)]./mean(Monocytes_Lucas_HCW),monocyte_sample_dayrand(where(end-49:end))'/p.M0];
figure
hold on 
violinplot(violin_mat_monocyte)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Monocytes')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


violin_mat_neutrophils = [Neutrophils_Lucas_moderate./mean(Neutrophils_Lucas_HCW), [neutrophil_sample_dayrand(where(1:50))'/p.N0;NaN(57-50,1)], [Neutrophils_Lucas_Severe;NaN(57-32,1)]./mean(Neutrophils_Lucas_HCW), [neutrophil_sample_dayrand(where(end-49:end))'/p.N0;NaN(57-50,1)]];
figure
hold on 
violinplot(violin_mat_neutrophils)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Neutrophils')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


violin_mat_IL6 = [10.^IL6_Lucas_moderate_individual'./mean(10.^IL6_Lucas_HCW_individual),[IL6_sample_dayrand(where(1:50))'/p.L_U_0;NaN(108-50,1)],[10.^IL6_Lucas_Severe_individual'./mean(10.^IL6_Lucas_HCW_individual);NaN(108-41,1)],[IL6_sample_dayrand(where(end-49:end))'/p.L_U_0;NaN(108-50,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('IL-6')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


%the maximum values for monocytes, neutrophils and IL6
violin_mat_monocyte = [[Monocytes_Lucas_moderate;NaN(2,1)], monocyte_sample_dayrand(where(1:50))'*1e9,[Monocytes_Lucas_Severe;NaN(50-29,1)],monocyte_sample_dayrand(where(end-49:end))'*1e9];
figure
hold on 
violinplot(violin_mat_monocyte)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Monocytes')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


violin_mat_neutrophils = [Neutrophils_Lucas_moderate, [neutrophil_sample_dayrand(where(1:50))'*1e9;NaN(57-50,1)], [Neutrophils_Lucas_Severe;NaN(57-32,1)], [neutrophil_sample_dayrand(where(end-49:end))'*1e9;NaN(57-50,1)]];
figure
hold on 
violinplot(violin_mat_neutrophils)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('Neutrophils')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')


violin_mat_IL6 = [10.^IL6_Lucas_moderate_individual',[IL6_sample_dayrand(where(1:50))';NaN(108-50,1)],[10.^IL6_Lucas_Severe_individual';NaN(108-41,1)],[IL6_sample_dayrand(where(end-49:end))';NaN(108-50,1)]];
figure
hold on 
violinplot(violin_mat_IL6)
ylabel('Change from baseline')
set(gca,'FontSize',18)
title('IL-6')
set(gca,'xticklabels',{'Mild (c)','Mild (p)','Severe (c)','Severe (p)'})
set(gca,'yscale','log')

STOP
%% Plot of IL-6 IFN and macs over time with shaded regions 

figure
hold on 
plot(time_grid,sol_IL6)
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('IL-6 (pg/ml)')

figure
hold on 
plot(time_grid,sol_IFN')
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('IFN (pg/ml)')

figure
hold on 
plot(time_grid,sol_macs_inflam')
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('Inflam macrophages (pg/ml)')

fig1 = figure
fig2 = figure
fig3 = figure

Psi_vec = linspace(min(Psi),max(Psi),10);

color_mat = [103,0,31
178,24,43;...
214,96,77;...
244,165,130;...
253,219,199;...
209,229,240;...
146,197,222;...
67,147,195;...
33,102,172;...
5,48,97]./255;%winter(10);

for ii = 1:200
   
    if Psi(ii) == min(Psi)
        Psi_color = 1;
    elseif Psi(ii) == max(Psi)
        Psi_color = 10;
    else
        Psi_color = find(Psi_vec>Psi(ii),1);
    end
    figure(fig1)
    hold on
   plot(time_grid,sol_IL6(ii,:),'Color',color_mat(Psi_color,:))
   
   
    figure(fig2)
    hold on
   plot(time_grid,sol_IFN(ii,:),'Color',color_mat(Psi_color,:))
   
   
    figure(fig3)
    hold on
   plot(time_grid,sol_macs_inflam(ii,:),'Color',color_mat(Psi_color,:))
    
end

figure(fig1)
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('IL-6 (pg/ml)')

figure(fig2)
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('IFN (pg/ml)')

figure(fig3)
xlabel('Time (days)')
xlim([0 21])
set(gca,'FontSize',16)
ylabel('Inflam macrophages (pg/ml)')
c=colorbar
colormap(color_mat)
set(c,'TickLabels',{'2','2.4','2.7','3','3.4','3.8'})

%%

%peaking = find(IFN_peak(where(51:end))<=4);

colmap = [213,62,79
244,109,67
253,174,97
254,224,139
255,255,191
230,245,152
171,221,164
102,194,165
50,136,189]/255;

IFN_peak_vec = linspace(1,10,10);
figure
hold on 
for ij = 1:length(patient)
   
    if IFN_peak(where(ij))==min(IFN_peak)
        colval = 2;
    else
        colval = find(IFN_peak_vec>=IFN_peak(where(ij)),1);
    end
    if colval==2
        colval
        IFN_peak(where(ij))
    end
    plot(ij,patient(where(ij),2),'.','Color',colmap(colval-1,:),'MarkerSize',30) 
end
xlabel('Patients')
ylabel('p_{M\Phi_I,L}')
set(gca,'FontSize',18)
[R,P] = corrcoef(IFN_peak(where),patient(where,2));
text(2,0.6,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)

colormap(colmap)
colorbar

figure
hold on 
plot(patient(where,3),'o') 
%plot(50+peaking,patient(where(50+peaking),3),'go')
xlabel('Patients')
ylabel('p_{L,M\Phi_I}')
set(gca,'FontSize',18)

figure
hold on 
plot(patient(where,4),'o') 
%plot(50+peaking,patient(where(50+peaking),4),'go')
xlabel('Patients')
ylabel('p_{F,I}')
set(gca,'FontSize',18)

figure
hold on 
plot(patient(where,7),'o') 
%plot(50+peaking,patient(where(50+peaking),7),'ko')
xlabel('Patients')
ylabel('p_{M,I}')
set(gca,'FontSize',18)

figure
hold on 
plot(patient(where,8),'o') 
%plot(50+peaking,patient(where(50+peaking),8),'go')
xlabel('Patients')
ylabel('\eta_{F,M\Phi}')
set(gca,'FontSize',18)


figure
hold on 
plot(patient(where,10),'o') 
%plot(50+peaking,patient(where(50+peaking),10),'go')
xlabel('Patients')
ylabel('\epsilon_{F,I}')
set(gca,'FontSize',18)


figure
hold on 
plot(patient(where,11),'o') 
%plot(50+peaking,patient(where(50+peaking),11),'go')
xlabel('Patients')
ylabel('p_{F,M}')
set(gca,'FontSize',18)

colmap = [247,252,253;...
224,236,244;...
191,211,230;...
158,188,218;...
140,150,198;...
140,107,177;...
136,65,157;...
129,15,124;...
77,0,75]/255;

IL6_grid = linspace(min(max_IL6),max(max_IL6),8);
figure
hold on 
for i = 1:length(patient)
    colindex = find(IL6_grid>=max_IL6(where(i)),1);
    plot3(sorted(i),max_T_cells(where(i))*1e9,max_IL6(where(i)),'.','Color',colmap(colindex,:),'LineWidth',14,'MarkerSize',40)
end

grid on 
colormap(colmap)
colorbar
set(gca,'FontSize',18)
xlabel('Inflam marker')
ylabel('Max T cells')
zlabel('Max IL-6')

colmap = [247,252,253;...
224,236,244;...
191,211,230;...
158,188,218;...
140,150,198;...
140,107,177;...
136,65,157;...
129,15,124;...
77,0,75]/255;

IL6_grid = linspace(min(max_IL6),max(max_IL6),8);
figure
hold on 
for i = 1:length(patient)
    colindex = find(IL6_grid>=max_IL6(where(i)),1);
    plot3(sorted(i),IFN_peak(where(i)),max_IL6(where(i)),'.','Color',colmap(colindex,:),'LineWidth',14,'MarkerSize',40)
end

grid on 
colormap(colmap)
colorbar
set(gca,'FontSize',18)
xlabel('Inflam marker')
ylabel('Max T cells')
zlabel('Max IL-6')


%% "correlation" matrices for IL-6 and tissue remaining

[sorted,where] = sort(real(minimum_tissue));
sol_tissue(jj,:) = deval(solstruc,time_grid,2)+deval(solstruc,time_grid,4);
sol_IL6(jj,:) = deval(solstruc,time_grid,11); 
    
IL6_color_grid = linspace(min(min(sol_IL6)),max(max(sol_IL6)),10);
tissue_color_grid = linspace(min(min(sol_tissue)),max(max(sol_tissue)),10);

col = [103,0,31
178,24,43;...
214,96,77;...
244,165,130;...
253,219,199;...
209,229,240;...
146,197,222;...
67,147,195;...
33,102,172;...
5,48,97]./255;%winter(10);

time_reduced = [find(time_grid>=1,1),find(time_grid>=2,1),find(time_grid>=3,1),...
    find(time_grid>=4,1),find(time_grid>=5,1),find(time_grid>=6,1),find(time_grid>=7,1),find(time_grid>=8,1),find(time_grid>=9,1),...
    find(time_grid>=10,1),find(time_grid>=11,1),find(time_grid>=12,1),find(time_grid>=13,1),find(time_grid>=14,1),find(time_grid>=15,1),...
    find(time_grid>=16,1),find(time_grid>=17,1),find(time_grid>=18,1),find(time_grid>=19,1),find(time_grid>=20,1)];
    
fig1 = figure
fig2 = figure
hold on
for ij = 1:length(patient)
    for jj = 1:20
        valIL6 = find(IL6_color_grid>=sol_IL6(where(ij),time_reduced(jj)),1);
        valtissue = find(tissue_color_grid>=sol_tissue(where(ij),time_reduced(jj)),1);
        
      figure(fig1)
      hold on
      fill([jj-1 jj-1 jj jj],[ij-1 ij ij ij-1],col(valIL6,:),'EdgeColor','none') 
       figure(fig2)
       hold on
       fill([jj-1 jj-1 jj jj],[ij-1 ij ij ij-1],col(valtissue,:),'EdgeColor','none') 
    end
end



%%

[sorted,where] = sort(minimum_tissue);

figure
waterfall_map = jet(10);
xgrid = linspace(1,200,200);

subplot(3,4,1)
hold on 
param_vec_sorted = patients(where,1);
param_grid = linspace(min(patients(:,1)),max(patients(:,1))*1.01,10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('\beta')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,2)
hold on 
param_vec_sorted = patients(where,2);
param_grid = linspace(min(patients(:,2)),max(patients(:,2)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('\epsilon_{F,I}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,3)
hold on 
param_vec_sorted = patients(where,3);
param_grid = linspace(min(patients(:,3)),max(patients(:,3)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('p_{F,M\Phi}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,4)
hold on 
param_vec_sorted = patients(where,4);
param_grid = linspace(min(patients(:,4)),max(patients(:,4)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('p_{F,M}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,5)
hold on 
param_vec_sorted = patients(where,5);
param_grid = linspace(min(patients(:,5)),max(patients(:,5)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('p_{F,I}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,6)
hold on 
param_vec_sorted = patients(where,6);
param_grid = linspace(min(patients(:,6)),max(patients(:,6)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('k_{B_F}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,7)
hold on 
param_vec_sorted = patients(where,7);
param_grid = linspace(min(patients(:,7)),max(patients(:,7)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('p_{M\Phi I,L}')   
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,8)
hold on 
param_vec_sorted = patients(where,8);
param_grid = linspace(min(patients(:,8)),max(patients(:,8)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('\eta_{F,I}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,9)
hold on 
param_vec_sorted = patients(where,9);
param_grid = linspace(min(patients(:,9)),max(patients(:,9)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('\epsilon_{L,T}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,10)
hold on 
param_vec_sorted = patients(where,10);
param_grid = linspace(min(patients(:,10)),max(patients(:,10)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('p_{M,I}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

subplot(3,4,11)
hold on 
param_vec_sorted = patients(where,11);
param_grid = linspace(min(patients(:,11)),max(patients(:,11)),10);
for i = 1:200
  
    param_val = param_vec_sorted(i);
    color_loc = find(param_grid>=param_val,1);
    
    bar(xgrid(i),sorted(i)./0.16*100,'FaceColor', waterfall_map(color_loc,:),'EdgeColor','none')

end
xlabel('Patient')
ylabel('Min %tissue remaining')
set(gca,'FontSize',12)
colormap(jet);
c = colorbar;
title('\eta_{F,M\Phi}')
ylim([min(sorted./0.16*100),max(sorted./0.16*100)])

%% Correlation matrix plot

colmap = [103,0,31
178,24,43;...
214,96,77;...
244,165,130;...
253,219,199;...
209,229,240;...
146,197,222;...
67,147,195;...
33,102,172;...
5,48,97]./255;%winter(10);
Rvalvec = linspace(0,1,10);
figure
hold on

patients = patient;
% remove_NAN_tissue = find(isnan(minimum_tissue)==1);
% minimum_tissue(remove_NAN_tissue) = [];
% patients(remove_NAN_tissue,:)=[];
% max_IL6(remove_NAN_tissue) = [];
% max_dead_tissue(remove_NAN_tissue) = [];
% IFN_exposure(remove_NAN_tissue) = [];
% IFN_peak(remove_NAN_tissue) = [];
% max_T_cells(remove_NAN_tissue) = [];
% max_infected_cells(remove_NAN_tissue) = [];
% max_inflam_macs(remove_NAN_tissue) = [];
% max_neutrophils(remove_NAN_tissue) = [];

for jk = 1:9
    R = corrcoef(patients(:,jk),minimum_tissue);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(1,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)

    R = corrcoef(patients(:,jk),max_IL6);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(1.5,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)

    R = corrcoef(patients(:,jk),max_dead_tissue);
    Rvalloc = find(Rvalvec>abs(R(1,2)),1);
    plot(2,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)

    R = corrcoef(patients(:,jk),IFN_exposure);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(2.5,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)
    
    R = corrcoef(patients(:,jk),IFN_peak);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(3,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)

    R = corrcoef(patients(:,jk),max_T_cells);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(3.5,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)
    
    R = corrcoef(patients(:,jk),max_infected_cells);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(4,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)
    
    R = corrcoef(patients(:,jk),max_inflam_macs);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(4.5,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)
    
    R = corrcoef(patients(:,jk),max_neutrophils);
    Rvalloc = find(Rvalvec>=abs(R(1,2)),1);
    plot(5,jk,'o','Color',colmap(Rvalloc,:),'LineWidth',13)
    
jk
end
colorbar
colormap(colmap)
%set(gca,'Colormap',flipud(colmap))
xlim([0.8 5])
set(gca,'Ytick',[1 2 3 4 5 6 7 8 9 10 11 12],'YtickLabels',{'p_{MPhi_I,L}','p_{L,MPhi}','p_{F,I}','p_{M,I}','\eta_{F,MPhi}','\epsilon_{F,I}','p_{M,F}'});%{'\beta','\epsilon_{F,I}','p_{F,MPhi}','p_{F,M}','p_{F,I}','k_{B_L}','p_{M\Phi I,L}','\eta_{F,I}','\epsilon_{L,T}','p_{M,I}','\eta_{F,M\Phi}','\tau_T'});
set(gca,'Xtick',[1 1.5 2 2.5 3 3.5 4 4.5 5],'XtickLabels',{'min(S+R)','max(L_U)','max(D)','Total IFN','IFN peak','max(T)','max(I)','max(MPhiI)','max(N)'});
    xtickangle(60)
    set(gca,'FontSize',13)
    ylim([0.5 9.5])