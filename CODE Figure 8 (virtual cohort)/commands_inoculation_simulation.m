% inoculation simulation


%load('patients_sim_annealing_severe6.mat')%2
load('patients_sim_annealing_severeNEW2.mat') % if no good, try 4 or no number
patients_full_cohort = patients_new;
patient = patients_new;
%load('patients_sim_annealing_MC.mat')
%patients_full_cohort = [patients_full_cohort;patient];

tspan = [0 21];
time_grid = linspace(tspan(1),tspan(2),1000);

p.p_F_IORIG = p.p_F_I;
p.eta_F_IORIG = p.eta_F_I;
p.p_M_IORIG = p.p_M_I;
p.eta_F_MPhiORIG = p.eta_F_MPhi;

p.I0 = 0;
dosage_vec = [0.1, 1, 4.5, 8, 12];

for kk = 1:5
    p.V0 = dosage_vec(kk);
    
    p.V0

    for jj = 1:200

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

        sol_virus(kk,jj,:) = deval(solstruc,time_grid,1);
        sol_tissue(kk,jj,:) = deval(solstruc,time_grid,2)+deval(solstruc,time_grid,4);
        sol_infected(kk,jj,:) = deval(solstruc,time_grid,3);
        sol_dead(kk,jj,:) = deval(solstruc,time_grid,5);
        sol_macs_res(kk,jj,:) = deval(solstruc,time_grid,6);
        sol_macs_inflam(kk,jj,:) = deval(solstruc,time_grid,7);
        sol_monocytes(kk,jj,:) = deval(solstruc,time_grid,8);
        sol_neutrophils(kk,jj,:) = deval(solstruc,time_grid,9);
        sol_Tcells(kk,jj,:) = deval(solstruc,time_grid,10);

        sol_IL6(kk,jj,:) = deval(solstruc,time_grid,11); 
        sol_IFN(kk,jj,:) = deval(solstruc,time_grid,17); 
        sol_GCSF(kk,jj,:) = deval(solstruc,time_grid,15);

        days14 = find(time>14,1);

        if isreal(sol)==1
            minimum_tissue(kk,jj) = min(sol(2,:)+sol(4,:));%+sol(3,:)
            max_IL6(kk,jj) = max(sol(11,:));
            max_bound_IL6(kk,jj) = max(sol(12,:));
            max_dead_tissue(kk,jj) = max(sol(5,:));
            max_inflam_macs(kk,jj) = max(sol(7,:));

            tgrid_vec = [];
            for i = 1:length(sol(17,:))-1
               tgrid_vec(i) = time(i+1)-time(i); 
            end
            IFN_exposure(kk,jj) = sum(sol(17,1:end-1).*tgrid_vec);

            peak_loc = find(sol(17,:)==max(sol(17,:)));
            IFN_peak(kk,jj) = time(peak_loc); %time of IFN peak
            max_T_cells(kk,jj) = max(sol(10,1:days14));    
            max_infected_cells(kk,jj) = max(sol(3,:));
            max_inflam_macs(kk,jj) = max(sol(7,:));
            max_neutrophils(kk,jj) = max(sol(9,:));
            max_monocytes(kk,jj) = max(sol(8,:));
        else
            minimum_tissue(kk,jj) = NaN;%+sol(3,:)
            max_IL6(kk,jj) = NaN;
            max_bound_IL6(kk,jj) = NaN;
            max_dead_tissue(kk,jj) = NaN;
            IFN_exposure(kk,jj) = NaN;
            IFN_peak(kk,jj) = NaN; %time of IFN peak
            max_T_cells(kk,jj) = NaN;    
            max_infected_cells(kk,jj) = NaN;
            max_inflam_macs(kk,jj) = NaN;
            max_neutrophils(kk,jj) = NaN;   
            max_monocytes(kk,jj) = NaN;
        end
        jj
    end
end

save('Patient responses_inoculation_NEW1.mat','time_grid','minimum_tissue','max_IL6','IFN_exposure','IFN_peak','max_T_cells','max_inflam_macs','max_neutrophils','max_monocytes','sol_virus','sol_tissue','sol_infected','sol_dead','sol_macs_res','sol_macs_inflam','sol_neutrophils','sol_monocytes','sol_Tcells','sol_IL6','sol_IFN','sol_GCSF')

STOP

%%
figure
hold on 

tspan=[0 21];
time_grid = linspace(tspan(1),tspan(2),1000);
p.Smax = 0.16;
[max_IL6ORIG max_neutrophilsORIG minimum_tissueORIG, max_inflam_macsORIG, IFN_exposureORIG, IFN_peakORIG, max_T_cellsORIG]=reevaluate_patients_special(time_grid,sol_virus,sol_tissue,sol_infected,sol_dead,sol_macs_res,sol_macs_inflam,sol_neutrophils, sol_monocytes,sol_Tcells, sol_IL6, sol_IFN, sol_GCSF,time_grid);

origPsi = (max_IL6ORIG)/(mean(max_IL6ORIG))+(max_neutrophilsORIG)/mean(max_neutrophilsORIG)+(p.Smax-minimum_tissueORIG)/mean(p.Smax-minimum_tissueORIG);

colvec = [237,248,177;...
199,233,180;...
127,205,187;...
65,182,196;...
29,145,192;...
34,94,168;...
37,52,148;...
8,29,88]/255;

Psi_vec = linspace(0,11,9);

for ii = 1:200
    Psi_v01 = (max_IL6(1,ii))/(nanmean(max_IL6(1,:)))+(max_neutrophils(1,ii))/nanmean(max_neutrophils(1,:))+(p.Smax-minimum_tissue(1,ii))/nanmean(p.Smax-minimum_tissue(1,:));
    Psi_v01
    if isnan(Psi_v01) == 0
        col = find(Psi_vec>=Psi_v01,1)-1;
        if isempty(col)==1
            col = 8;
        end
        plot(origPsi(ii),max_IL6(1,ii),'o','Color',colvec(col,:),'LineWidth',2)
        Psi_v01_vec(ii) = Psi_v01;
    else
        Psi_v01_vec(ii) = NaN;
    end
end
xlabel('Original Psi')
ylabel('Max IL-6')
colormap(colvec)
c = colorbar;
title('V0 = 0.1')
set(gca,'FontSize',18)

for ii = 1:200
    Psi_v02 = (max_IL6(2,ii))/(nanmean(max_IL6(2,:)))+(max_neutrophils(2,ii))/nanmean(max_neutrophils(2,:))+(p.Smax-minimum_tissue(2,ii))/nanmean(p.Smax-minimum_tissue(2,:));
    Psi_v02
    if isnan(Psi_v02) == 0
        col = find(Psi_vec>=Psi_v02,1)-1;
        if isempty(col)==1
            col = 8;
        end
        plot(origPsi(ii),max_IL6(2,ii),'o','Color',colvec(col,:),'LineWidth',2)        
        Psi_v02_vec(ii) = Psi_v02;
    else
        Psi_v02_vec(ii) = NaN;
    end
end
xlabel('Original Psi')
ylabel('Max IL-6')
colormap(colvec)
c = colorbar;
title('V0 = 1')

for ii = 1:200
    Psi_v03 = (max_IL6(3,ii))/(nanmean(max_IL6(3,:)))+(max_neutrophils(3,ii))/nanmean(max_neutrophils(3,:))+(p.Smax-minimum_tissue(3,ii))/nanmean(p.Smax-minimum_tissue(3,:));
    Psi_v03
    if isnan(Psi_v03) == 0
        col = find(Psi_vec>=Psi_v03,1)-1;
        if isempty(col)==1
            col = 8;
        end
        plot(origPsi(ii),max_IL6(3,ii),'o','Color',colvec(col,:),'LineWidth',2)
        Psi_v03_vec(ii) = Psi_v03;
    else
        Psi_v03_vec(ii) = NaN;
    end
end
xlabel('Original Psi')
ylabel('Max IL-6')
colormap(colvec)
c = colorbar;
title('V0 = 8')

for ii = 1:200
    Psi_v04 = (max_IL6(4,ii))/(nanmean(max_IL6(4,:)))+(max_neutrophils(4,ii))/nanmean(max_neutrophils(4,:))+(p.Smax-minimum_tissue(4,ii))/nanmean(p.Smax-minimum_tissue(4,:));
    Psi_v04
    if isnan(Psi_v04) == 0
        col = find(Psi_vec>=Psi_v04,1)-1;
        if isempty(col)==1
            col = 8;
        end
        plot(origPsi(ii),max_IL6(4,ii),'o','Color',colvec(col,:),'LineWidth',2)
        Psi_v04_vec(ii) = Psi_v04;
    else
        Psi_v04_vec(ii) = NaN;
    end
end
xlabel('Original Psi')
ylabel('Max IL-6')
colormap(colvec)
c = colorbar;
title('V0 = 12')
set(c,'ticks',linspace(0,1,9),'ticklabels',{'0','1.4','2.8','4.1','5.5','6.9','8.3','9.6','11'})


for ii = 1:200
    Psi_v05 = (max_IL6(5,ii))/(nanmean(max_IL6(5,:)))+(max_neutrophils(5,ii))/nanmean(max_neutrophils(5,:))+(p.Smax-minimum_tissue(5,ii))/nanmean(p.Smax-minimum_tissue(5,:));
    Psi_v05
    if isnan(Psi_v05) == 0
        col = find(Psi_vec>=Psi_v05,1)-1;
        if isempty(col)==1
            col = 8;
        end
        plot(origPsi(ii),max_IL6(5,ii),'o','Color',colvec(col,:),'LineWidth',2)
        Psi_v05_vec(ii) = Psi_v05;
    else
        Psi_v05_vec(ii) = NaN;
    end
end
xlabel('Original Psi')
ylabel('Max IL-6')
colormap(colvec)
c = colorbar;
title('V0 = 12')
set(c,'ticks',linspace(0,1,9),'ticklabels',{'0','1.4','2.8','4.1','5.5','6.9','8.3','9.6','11'})

%%

figure
hold on 
plot(origPsi,Psi_v01_vec,'o','Color',[213,62,79]/255)
plot(origPsi,Psi_v02_vec,'o','Color',[171,221,164]/255)
plot(origPsi,Psi_v04_vec,'o','Color',[50,136,189]/255)
l1 = plot(0,0,'o','Color',[213,62,79]/255,'LineWidth',1.5)
l2 = plot(0,0,'o','Color',[171,221,164]/255,'LineWidth',1.5)
l3 = plot(0,0,'o','Color',[50,136,189]/255,'LineWidth',1.5)
xlabel('Inflammation marker \Psi at V_0 = 4.5')
ylabel({'Inflammation marker \Psi at', 'V_0 = 0.1, 1 and 8'})
legend([l1 l2 l3],{'V_0 = 0.1','V_0 = 1','V_0 = 8'})
set(gca,'FontSize',18)
xlim([2 4])
ylim([2 4])
nonan_Psi_v01_Vec = Psi_v01_vec;
nonan_Psi_v01_Vec(isnan(Psi_v01_vec))=[];
nonan_origPsi = origPsi;
nonan_origPsi(isnan(Psi_v01_vec))=[];
[R,P] = corrcoef(nonan_origPsi,nonan_Psi_v01_Vec);
text(2.05,2.2,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
[R,P] = corrcoef(nonan_origPsi,nonan_Psi_v02_Vec);
text(2.05,2.2,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)
[R,P] = corrcoef(nonan_origPsi,nonan_Psi_v03_Vec);
text(2.05,2.2,['R = ' num2str(R(1,2)) ', p = ' num2str(P(1,2))],'FontSize',14)

figure
hold on 
plot(Psi_v03_vec,Psi_v01_vec,'o')
plot(Psi_v03_vec,Psi_v02_vec,'o')
%plot(origPsi,Psi_v03_vec,'o:')
plot(Psi_v03_vec,Psi_v04_vec,'o')
xlabel('Inflammation marker \Psi at V_0 = 4.5')
ylabel('Inflammation marker \Psi')
legend('V_0 = 0.1','V_0 = 1','V_0 = 8')
set(gca,'FontSize',18)

figure
hold on 
plot(origPsi(1:50),max_IL6(1,1:50),'o:')
plot(origPsi(1:50),max_IL6(2,1:50),'o:')
plot(origPsi(1:50),max_IL6(3,1:50),'o:')
plot(origPsi(1:50),max_IL6(4,1:50),'o:')
xlabel('Original Psi')
ylabel('Max IL-6( pg/ml)')
legend('V0 = 0.1','V0 = 1','V0 = 8','V0 = 12')


figure
hold on 
plot(origPsi(1:50),max_neutrophils(1,1:50),'o','LineWidth',2)
plot(origPsi(1:50),max_neutrophils(2,1:50),'o','LineWidth',2)
plot(origPsi(1:50),max_neutrophils(3,1:50),'o','LineWidth',2)
plot(origPsi(1:50),max_neutrophils(4,1:50),'o','LineWidth',2)
xlabel('Original Psi')
ylabel('Max neutrophils( cells/ml)')
legend('V0 = 0.1','V0 = 1','V0 = 8','V0 = 12')
%%

colvec = [208,209,230
166,189,219
103,169,207
54,144,192
2,129,138
1,108,89
1,70,54]/255;

neut_vec = linspace(min(max_neutrophilsORIG),max(max_neutrophilsORIG),8);

psi_vecO = linspace(min(origPsi),max(origPsi),8);

figure
hold on 
for ii = 1:length(patient)
    col = find(psi_vecO>origPsi(ii),1)-1;%neut_vec>max_neutrophilsORIG(ii),1)-1;
    if isempty(col) == 1
       col = 7; 
    end
   % plot([1 2 3 4 5],max_neutrophils(:,ii),':','Color',colvec(col,:),'LineWidth',0.5)
    plot([1 2 3 4],max_neutrophils(:,ii),'o','Color',colvec(col,:),'LineWidth',1.5)
end
xlabel('Exposure dose (log_{10}(copies/ml))')
ylabel('Max neutrophils (cells/ml)')
set(gca,'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original neuts')
colormap(colvec)
set(gca,'FontSize',18)
c = colorbar
set(c,'ticks',linspace(0,1,8),'ticklabels',{'5.3','5.8','6.3','6.7','7.2','7.6','8.0','8.5'}) % FIX
saveas(gcf,'Fig_8H.fig');
saveas(gcf,'Fig_8H.png');

il6_vec = linspace(min(max_IL6ORIG),max(max_IL6ORIG),8);

figure
hold on 
for ii = 1:length(patient)
        col = find(psi_vecO>origPsi(ii),1)-1;%il6_vec>max_IL6ORIG(ii),1)-1;
    if isempty(col) == 1
       col = 7; 
    end
   % plot([1 2 3 4 5],max_IL6(:,ii),':','Color',colvec(col,:),'LineWidth',0.5)
    plot([1 2 3 4],max_IL6(:,ii),'o','Color',colvec(col,:),'LineWidth',1.5)
end
xlabel('Exposure dose (log_{10}(copies/ml))')
ylabel('Max IL-6 (pg/ml)')
set(gca,'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original max IL6')
colormap(colvec)
set(gca,'FontSize',18)
ylim([18 81])
c = colorbar
set(c,'ticks',linspace(0,1,8),'ticklabels',{'21','30','40','48','57','67','76','85'}) % FIX
saveas(gcf,'Fig_8HH.fig');
saveas(gcf,'Fig_8HH.png');

psi_vec = linspace(min(origPsi),max(origPsi),8);

figure
hold on 
for ii = 1:length(patient)
        col = find(psi_vecO>origPsi(ii),1)-1;%psi_vec>origPsi(ii),1)-1;
    if isempty(col) == 1
       col = 7; 
    end
    %plot([1 2 3 4 5],[Psi_v01_vec(ii),Psi_v02_vec(ii),Psi_v03_vec(ii),Psi_v04_vec(ii),Psi_v05_vec(ii)],':','Color',colvec(col,:),'LineWidth',0.5)
    plot([1 2 3 4],[Psi_v01_vec(ii),Psi_v02_vec(ii),Psi_v03_vec(ii),Psi_v04_vec(ii)],'o','Color',colvec(col,:),'LineWidth',1.5)   %,Psi_v05_vec(ii)
end
xlabel('Dosage (log_{10}(copies/ml))')
ylabel('Inflammation marker \Psi')
set(gca,'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original Psi')
colormap(colvec)
c = colorbar
set(gca,'FontSize',18)
set(c,'ticks',linspace(0,1,8),'ticklabels',{'2.2','2.4','2.6','2.9','3','3.3','3.6','3.7'}) % FIX
saveas(gcf,'Fig_8HHH.fig');
saveas(gcf,'Fig_8HHH.png');

%%


colvec = [208,209,230
166,189,219
103,169,207
54,144,192
2,129,138
1,108,89
1,70,54]/255;

neut_vec = linspace(min(min(max_neutrophils)),max(max(max_neutrophils)),8);
[sorted,where] = sort(origPsi);


figure
hold on 
for ii = 1:length(patient)
    ID = where(ii);
    col1 = find(neut_vec>max_neutrophils(1,ID),1)-1;
    col2 = find(neut_vec>max_neutrophils(2,ID),1)-1;
    col3 = find(neut_vec>max_neutrophils(3,ID),1)-1;
    col4 = find(neut_vec>max_neutrophils(4,ID),1)-1;
    col5 = find(neut_vec>max_neutrophils(5,ID),1)-1;
    if isempty(col1) == 1
       col1 = 7; 
    elseif isempty(col2) == 1
        col2 = 7;
    elseif isempty(col3) == 1
        col3 = 7;
    elseif isempty(col4) == 1
        col4 = 7;
    elseif isempty(col5) == 1
        col5 = 7;
    end
    fill([0 0 1 1],[ii-1 ii ii ii-1],colvec(col1,:),'EdgeColor','none') 
    fill([1 1 2 2],[ii-1 ii ii ii-1],colvec(col2,:),'EdgeColor','none') 
    fill([2 2 3 3],[ii-1 ii ii ii-1],colvec(col3,:),'EdgeColor','none') 
    fill([3 3 4 4],[ii-1 ii ii ii-1],colvec(col4,:),'EdgeColor','none') 
    fill([4 4 5 5],[ii-1 ii ii ii-1],colvec(col5,:),'EdgeColor','none') 
end
hold on 

xlabel('Dosage (log_{10}(copies/ml))')
ylabel('Patients (ordered by \Psi)')
set(gca,'xtick',[0.5 1.5 2.5 3.5 4.5],'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original neuts')
colormap(colvec)
set(gca,'FontSize',18)
c = colorbar
set(c,'ticks',linspace(0,1,8),'ticklabels',{'5.3','5.8','6.3','6.7','7.2','7.6','8.0','8.5'}) % FIX
xlim([1 5])
saveas(gcf,'Fig_8H_1.fig');
saveas(gcf,'Fig_8H_1.png');
 
il6_vec = linspace(min(min(max_IL6)),max(max(max_IL6)),8);
figure
hold on 
for ii = 1:length(patient)
    ID = where(ii);
    col1 = find(il6_vec>max_IL6(1,ID),1)-1;
    col2 = find(il6_vec>max_IL6(2,ID),1)-1;
    col3 = find(il6_vec>max_IL6(3,ID),1)-1;
    col4 = find(il6_vec>max_IL6(4,ID),1)-1;
    col5 = find(il6_vec>max_IL6(5,ID),1)-1;
    if isempty(col1) == 1
       col1 = 7; 
    elseif isempty(col2) == 1
        col2 = 7;
    elseif isempty(col3) == 1
        col3 = 7;
    elseif isempty(col4) == 1
        col4 = 7;
    elseif isempty(col5) == 1
        col5 = 7;
    end
    fill([0 0 1 1],[ii-1 ii ii ii-1],colvec(col1,:),'EdgeColor','none') 
    fill([1 1 2 2],[ii-1 ii ii ii-1],colvec(col2,:),'EdgeColor','none') 
    fill([2 2 3 3],[ii-1 ii ii ii-1],colvec(col3,:),'EdgeColor','none') 
    fill([3 3 4 4],[ii-1 ii ii ii-1],colvec(col4,:),'EdgeColor','none') 
    fill([4 4 5 5],[ii-1 ii ii ii-1],colvec(col5,:),'EdgeColor','none') 
end
hold on 
xlabel('Dosage (log_{10}(copies/ml))')
ylabel('Patients (ordered by \Psi)')
set(gca,'xtick',[0.5 1.5 2.5 3.5 4.5],'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original neuts')
colormap(colvec)
set(gca,'FontSize',18)
c = colorbar
set(c,'ticks',linspace(0,1,8),'ticklabels',{'21','30','40','48','57','67','76','85'}) % FIX
xlim([1 5])
saveas(gcf,'Fig_8HH_1.fig');
saveas(gcf,'Fig_8HH_1.png');

psi_mat = [Psi_v01_vec;Psi_v02_vec;Psi_v03_vec;Psi_v04_vec;Psi_v05_vec];
psi_vec = linspace(min(min([Psi_v01_vec;Psi_v02_vec;Psi_v03_vec;Psi_v04_vec;Psi_v05_vec])),max(max([Psi_v01_vec;Psi_v02_vec;Psi_v03_vec;Psi_v04_vec;Psi_v05_vec])),8);

figure
hold on 
for ii = 1:length(patient)
    ID = where(ii);
    col1 = find(psi_vec>psi_mat(1,ID),1)-1;
    col2 = find(psi_vec>psi_mat(2,ID),1)-1;
    col3 = find(psi_vec>psi_mat(3,ID),1)-1;
    col4 = find(psi_vec>psi_mat(4,ID),1)-1;
    col5 = find(psi_vec>psi_mat(5,ID),1)-1;
    if isempty(col1) == 1
       col1 = 7; 
    elseif isempty(col2) == 1
        col2 = 7;
    elseif isempty(col3) == 1
        col3 = 7;
    elseif isempty(col4) == 1
        col4 = 7;
    elseif isempty(col5) == 1
        col5 = 7;
    end
    fill([0 0 1 1],[ii-1 ii ii ii-1],colvec(col1,:),'EdgeColor','none') 
    fill([1 1 2 2],[ii-1 ii ii ii-1],colvec(col2,:),'EdgeColor','none') 
    fill([2 2 3 3],[ii-1 ii ii ii-1],colvec(col3,:),'EdgeColor','none') 
    fill([3 3 4 4],[ii-1 ii ii ii-1],colvec(col4,:),'EdgeColor','none') 
    fill([4 4 5 5],[ii-1 ii ii ii-1],colvec(col5,:),'EdgeColor','none') 
end
hold on 
xlabel('Dosage (log_{10}(copies/ml))')
ylabel('Patients (ordered by \Psi)')
set(gca,'xtick',[0.5 1.5 2.5 3.5 4.5],'xticklabels',{'0.1','1','4.5','8','12'})
%title('colorbar for original neuts')
colormap(colvec)
set(gca,'FontSize',18)
c = colorbar
set(c,'ticks',linspace(0,1,8),'ticklabels',{'2.2','2.4','2.6','2.9','3','3.3','3.6','3.7'}) % FIX
xlim([1 5])

saveas(gcf,'Fig_8HHH_1.fig');
saveas(gcf,'Fig_8HHH_1.png');