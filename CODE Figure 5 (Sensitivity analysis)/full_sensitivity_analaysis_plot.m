function full_sensitivity_analaysis_plot(mega_metric_matrix,mega_metric_matrix2)

fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
fig5 = figure;
fig6 = figure;
fig7 = figure;
fig8 = figure;

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

col_days = [64,0,75;...
118,42,131;...
153,112,171;...
194,165,207;...
231,212,232;...
247,247,247;...
217,240,211;...
166,219,160;...
90,174,97;...
27,120,55;...
0,68,27]/255;

comp_mat = [mega_metric_matrix;mega_metric_matrix2];

% 
figure(fig1)
hold on

for j = 1:76
    
    jval = j;
    
    if rem(j,2)==1%j is odd
        figure(fig1)
        hold on
        plotting_heatmap(mega_metric_matrix,1,8,jval,col,comp_mat)
        figure(fig3)
        hold on
        plotting_heatmap(mega_metric_matrix,9,10,jval,col_days,comp_mat)
    else % j is even
        figure(fig2)
        hold on
        plotting_heatmap(mega_metric_matrix,1,8,jval,col,comp_mat)
        figure(fig4)
        hold on
        plotting_heatmap(mega_metric_matrix,9,10,jval,col_days,comp_mat)
    end

end

%---------------------------------------------
figure(fig1)
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''}) %positive
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
xtickangle(45)
set(gca,'FontSize',14)
ylim([0 38])
xlim([0 8])

figure(fig2)
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
set(gca,'Ytick',linspace(1,38,38),'Yticklabel',{'p','\epsilon_{F,I}','d_V',...
    '\delta_{V,M\Phi}','\delta_{V,N}','\beta','\lambda_S','\delta_N',...
   'IC_{50,N}','d_I','\delta_{I,M\Phi}','\delta_{I,T}','\delta_{M\Phi,D}','\delta_{D,M\Phi}',...
   'a_{I,M\Phi}','\epsilon_{V,M\Phi}','\epsilon_{G,M\Phi}','p_{M\Phi,G}','\epsilon_{L,M\Phi}',...
   'p_{M\Phi,L}','\epsilon_{G,M}',...
    '\psi_M^{max}','p_{M,I}','\epsilon_{I,M}','\psi_N^{max}','\epsilon_{C,N}',...
    'p_{N,L}','\epsilon_{L,N}','\epsilon_{L,T}',...
    'p_{T,I}','\epsilon_{F,T}','p_{T,F}','V_0','S_0','M\Phi_{R,0}','M_0','N_0','T_0'}) % negative
xtickangle(45)
set(gca,'FontSize',14)
ylim([0.5 38.5])
xlim([0 8])


figure(fig3)
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''})
set(gca,'Xtick',[8.5,8.8,9.5],'Xticklabel',{'Time tissue','under 30%','Peak IFN time'})%
xtickangle(45)
set(gca,'FontSize',14)
ylim([0 38])

figure(fig4)
set(gca,'Xtick',[8.5,8.8,9.5],'Xticklabel',{'Time tissue','under 30%','Peak IFN time'})%
set(gca,'Ytick',linspace(1,38,38),'Yticklabel',{'p','\epsilon_{F,I}','d_V',...
    '\delta_{V,M\Phi}','\delta_{V,N}','\beta','\lambda_S','\delta_N',...
   'IC_{50,N}','d_I','\delta_{I,M\Phi}','\delta_{I,T}','\delta_{M\Phi,D}','\delta_{D,M\Phi}',...
   'a_{I,M\Phi}','\epsilon_{V,M\Phi}','\epsilon_{G,M\Phi}','p_{M\Phi,G}','\epsilon_{L,M\Phi}',...
   'p_{M\Phi,L}','\epsilon_{G,M}',...
    '\psi_M^{max}','p_{M,I}','\epsilon_{I,M}','\psi_N^{max}','\epsilon_{C,N}',...
    'p_{N,L}','\epsilon_{L,N}','\epsilon_{L,T}',...
    'p_{T,I}','\epsilon_{F,T}','p_{T,F}','V_0','S_0','M\Phi_{R,0}','M_0','N_0','T_0'})
    
xtickangle(45)
set(gca,'FontSize',14)
ylim([0.5 38.5])


%-----------------------------------------------------

for j = 1:76
    
    jval = j;
    
    if rem(j,2)==1%j is odd
        figure(fig5)
        hold on
        plotting_heatmap(mega_metric_matrix2,1,8,jval,col,comp_mat)
        figure(fig7)
        hold on
        plotting_heatmap(mega_metric_matrix2,9,10,jval,col_days,comp_mat)
    else % j is even
        figure(fig6)
        hold on
        plotting_heatmap(mega_metric_matrix2,1,8,jval,col,comp_mat)
        figure(fig8)
        hold on
        plotting_heatmap(mega_metric_matrix2,9,10,jval,col_days,comp_mat)
    end

end
%---------------------------------------------
figure(fig5)
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''}) %positive
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
xtickangle(45)
set(gca,'FontSize',14)
ylim([0 38])
xlim([0 8])

figure(fig6)
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
set(gca,'Ytick',linspace(1,38,38),'Yticklabel',{'p_{L,I}',...
    '\eta_{L,I}','p_{L,_M\Phi}','\eta_{L,M\Phi}','p_{L,_M}','\eta_{L,M}',...
    'k_{lin}_L','k_B_L','k_U_L','k_{int}_L','p_{G,M\Phi}','\eta_{G,M\Phi}',...
    'p_{G,M}','\eta_{G,M}','k_{lin}_G','k_B_G','k_U_G','k_{int}_G','p_{C,M}','\eta_{C,M}','k_{lin}_C','k_B_C','k_U_C','k_{int}_C',...
    'p_{F,I}','\eta_{F,I}','p_{F,MPhi}','\eta_{F,MPhi}',...
    'p_{F,M}','\eta_{F,M}','k_{lin}_F','k_B_F','k_U_F','k_{int}_F',...
   'L_{U,0}','G_{U,0}','C_{U,0}','F_{U,0}'}) % negative
xtickangle(45)
set(gca,'FontSize',14)
ylim([0.5 38.5])
xlim([0 8])


figure(fig7)
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''})
set(gca,'Xtick',[8.5,8.8,9.5],'Xticklabel',{'Time tissue','under 30%','Peak IFN time'})%
xtickangle(45)
set(gca,'FontSize',14)
ylim([0 38])

figure(fig8)
set(gca,'Xtick',[8.5,8.8,9.5],'Xticklabel',{'Time tissue','under 30%','Peak IFN time'})%
set(gca,'Ytick',linspace(1,38,38),'Yticklabel',{'p_{L,I}',...
    '\eta_{L,I}','p_{L,_M\Phi}','\eta_{L,M\Phi}','p_{L,_M}','\eta_{L,M}',...
    'k_{lin}_L','k_B_L','k_U_L','k_{int}_L','p_{G,M\Phi}','\eta_{G,M\Phi}',...
    'p_{G,M}','\eta_{G,M}','k_{lin}_G','k_B_G','k_U_G','k_{int}_G','p_{C,M}','\eta_{C,M}','k_{lin}_C','k_B_C','k_U_C','k_{int}_C',...
    'p_{F,I}','\eta_{F,I}','p_{F,MPhi}','\eta_{F,MPhi}',...
    'p_{F,M}','\eta_{F,M}','k_{lin}_F','k_B_F','k_U_F','k_{int}_F',...
   'L_{U,0}','G_{U,0}','C_{U,0}','F_{U,0}'}) % negative
xtickangle(45)
set(gca,'FontSize',14)
ylim([0.5 38.5])

end
