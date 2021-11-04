%%
tspan = [0 30];

[timeM,solM] = COVID_IMMUNE_MODEL(p,tspan);

%full_model_plotter(timeM,solM,p)
% 
% %Parameters to delay IFN/get severe dynamics
% p.p_F_I = 0.002;%0.0002
% %p.eta_F_I = 1;%16
% p.eta_F_MPhi = 1e-4*2;%1e-3
% p.p_M_I = 2*0.6;

p.a_I_MPhi = 1100*0.5;
p.p_MPhi_I_G = 1.68*0.5;
p.p_MPhi_I_L = 1.68*0.5;
p = Homeostasis_calculations(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end

[timeR,solR] = COVID_IMMUNE_MODEL(p,tspan);

p.MPhi_I_0 = 0;

p.a_I_MPhi = 0;
p.p_MPhi_I_G = 0;
p.p_MPhi_I_L = 0;

[timeS,solS] = COVID_IMMUNE_MODEL_macs(p,tspan);

timeMACS_deval = linspace(tspan(1),tspan(2),1e3);
solMMACS_deval = interp1(timeM,solM',timeMACS_deval);
solSMACS_deval = interp1(timeS,solS',timeMACS_deval);
solRMACS_deval = interp1(timeR,solR',timeMACS_deval);

mild_neut_changes_plotting(timeMACS_deval,solMMACS_deval',timeMACS_deval,solSMACS_deval',timeMACS_deval,solRMACS_deval',p)

p.a_I_MPhi = 1100;
p.p_MPhi_I_G = 1.68;
p.p_MPhi_I_L = 1.68;

%%---------------------------
tspan = [0 30];

[timeM,solM] = COVID_IMMUNE_MODEL(p,tspan);

%full_model_plotter(timeM,solM,p)
% 
% %Parameters to delay IFN/get severe dynamics
% p.p_F_I = 0.002;%0.0002
% %p.eta_F_I = 1;%16
% p.eta_F_MPhi = 1e-4*2;%1e-3
% p.p_M_I = 2*0.6;

p.p_M_I = 0.22*0.5;
p.psi_M_max = 0.133*0.5;% p.N_prod_star = 0.21     11.5451
p.M0 = 1e-4*0.5;    %4.000000000000000e-04
p = Homeostasis_calculations(p);   

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end

[timeR,solR] = COVID_IMMUNE_MODEL(p,tspan);

p.M0 = 0;

p.p_M_I = 0;
p.M_prod_star = 0;
p.psi_M_max = 0;
%p.C_U_0 = 0;
%p.C_B_0 = 0;

p.A_L = p.stoch*p.MM_L/p.avo*(p.R_L_T)*(1/5000)*10^9*1e12;% IL-6

%p = Homeostasis_calculations_mono(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end
%----------------------------------------------------------------------
p.M0 = 0;
p.p_M_I = 0;
p.M_prod_star = 0;
p.psi_M_max = 0;
p.d_V_spec = 2;

[timeS,solS] = COVID_IMMUNE_MODEL_mono(p,tspan);

timeMONOS_deval = linspace(tspan(1),tspan(2),1e3);
solMMONOS_deval = interp1(timeM,solM',timeMONOS_deval);
solSMONOS_deval = interp1(timeS,solS',timeMONOS_deval);
solRMONOS_deval = interp1(timeR,solR',timeMONOS_deval);

mild_neut_changes_plotting(timeMONOS_deval,solMMONOS_deval',timeMONOS_deval,solSMONOS_deval',timeMONOS_deval,solRMONOS_deval',p)

p.p_M_I = 0.22;
p.psi_M_max = 11.5451;
p.M0 = 4.000000000000000e-04;
p.M_prod_star = 0.133489926337691;
p.d_V_spec = 0;


%%


tspan = [0 30];

[timeM,solM] = COVID_IMMUNE_MODEL(p,tspan);

%full_model_plotter(timeM,solM,p)
% 
% %Parameters to delay IFN/get severe dynamics
% p.p_F_I = 0.002;%0.0002
% %p.eta_F_I = 1;%16
% p.eta_F_MPhi = 1e-4*2;%1e-3
% p.p_M_I = 2*0.6;

p.psi_N_max = 0.5;% p.N_prod_star = 0.21     4.13350
p.p_N_L = 1e-5;  %    0.212883162468174
p.N0 = 0.005*0.5;  %  0.00526

p = Homeostasis_calculations(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end

[timeR,solR] = COVID_IMMUNE_MODEL(p,tspan);

p.N0 = 0;
p.N_prod_star = 0;
p.psi_N_max = 0;
p.p_N_L = 0;
p.C_U_0 = 0;
p.C_B_0 = 0;

p.A_L = p.stoch*p.MM_L/p.avo*(p.R_L_M+p.R_L_T)*(1/5000)*10^9*1e12;% IL-6
p.A_G = p.stoch*p.MM_G/p.avo*(p.R_G_M)*(1/5000)*10^9*1e12;% GM-CSF
p.A_C = 0;%2*p.MM_C/p.avo*p.R_C_N*(1/5000)*10^9*1e9;% G-CSF
p.A_F = p.stoch*p.MM_F/p.avo*(p.R_F_T+p.R_F_I)*(1/5000)*10^9*1e12;% IFN

p = Homeostasis_calculations_neuts(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end
%----------------------------------------------------------------------
p.N0 = 0;
p.N_prod_star = 0;
p.psi_N_max = 0;
p.p_N_L = 0;


p.d_V_spec = 2;


[timeS,solS] = COVID_IMMUNE_MODEL_neuts(p,tspan);

timeNEUTS_deval = linspace(tspan(1),tspan(2),1e3);
solMNEUTS_deval = interp1(timeM,solM',timeNEUTS_deval);
solSNEUTS_deval = interp1(timeS,solS',timeNEUTS_deval);
solRNEUTS_deval = interp1(timeR,solR',timeNEUTS_deval);

mild_neut_changes_plotting(timeNEUTS_deval,solMNEUTS_deval',timeNEUTS_deval,solSNEUTS_deval',timeNEUTS_deval,solRNEUTS_deval',p)


p.psi_N_max = 4.13350;
p.p_N_L = 0.212883162468174;
p.N0 =0.00526;
p.N_prod_star = 0.212883162468174;
p.C_U_0 = 0.025;
p.C_B_0 = 6.503340242271287e-10;
p.d_V_spec = 0;


plot_monos_neuts_macs(timeNEUTS_deval,solMMACS_deval',solSMACS_deval',solRMACS_deval',solMMONOS_deval',solSMONOS_deval',solRMONOS_deval',solMNEUTS_deval',solSNEUTS_deval',solRNEUTS_deval',p)



