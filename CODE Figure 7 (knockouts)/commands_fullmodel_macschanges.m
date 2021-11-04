%simulating full model

format long

p = load_parameters();

% % ------- new changes to il-6 effector concentration 20th Aug
% p.eps_L_N = 5.7197*1e-5;
% p.eps_G_M = 3.950758130451212e-07;
% % -------

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
%p.p_G_M = 7.7*1e5;
%p.eps_G_MPhi = 1.846806147120571e-05;%p.eps_G_MPhi/(p.G_U_0/p.G_B_0)

%modulated IFN parameters to get desired IFN dynamics in the IFN submodel
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;

p.p_F_M = 3.5600*1e4*1e-4; %changed

%fixing IL-6 range
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;%1e-2*1.3;%1102.9/1e7;% HERE!!!!!-
p.eps_G_MPhi = 2664.5/1e5;%1e-3*1.3;%2664.5/1e7;%%1.846806147120571e-05;

p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;%0.5*1e-3/0.5;

p.p_T_I = 0.01*0.8; %changed from 1
p.del_I_T = 238;%50;%50;238 got from T cell killing paper

p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
%
%p.tau_T = 6;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; % CHANGED FROM 5
p.eps_I_M = 0.11;
p.del_V_N = 768*3;%*1.5;
p.del_V_MPhi = 768*100;%10000*4*10;%*1.5;

%p.d_MPhi_I = 1;
%p.k_lin_C = 16;%

%p.d_D = 21;
%p.d_M = p.d_N;
%p.p_MPhi_I_G = 0.42*2;
%p.p_MPhi_I_L = 0.42*2;

p.p_G_M = 1.234*1e3*1e-1;
%p.p_G_MPhi_I = 1232.4;

p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%p.eta_G_MPhi = 1e-3;

% XXXXXXXXXXXXXXX
% p.d_V = 18.94*1.2;
% p.d_V_spec = 8;
% p.d_I = 0.0144*6;
% p.phat = 741.2*2.5;
% p.beta = 0.289*0.6;
% p.del_V_N = 768*3*2;%*1.5;
% p.del_V_MPhi = 768*100*2;%10000*4*10;%*1.5;
% XXXXXXXXXXXXXXX

p = Homeostasis_calculations(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end
%----------------------------------------------------------------------

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

time_deval = linspace(tspan(1),tspan(2),1e3);
solM_deval = interp1(timeM,solM',time_deval);
solS_deval = interp1(timeS,solS',time_deval);
solR_deval = interp1(timeR,solR',time_deval);

mild_neut_changes_plotting(time_deval,solM_deval',time_deval,solS_deval',time_deval,solR_deval',p)
