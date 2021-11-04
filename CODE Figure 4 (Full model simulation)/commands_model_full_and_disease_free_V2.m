%simulating full model

format long

p = load_parameters();

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;
p.p_F_M = 3.5600*1e4*1e-4; 
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;
p.eps_G_MPhi = 2664.5/1e5;
p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;
p.p_T_I = 0.01*0.8;
p.del_I_T = 238;
p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; 
p.eps_I_M = 0.11;
p.del_V_N = 768*3;
p.del_V_MPhi = 768*100;
p.p_G_M = 1.234*1e3*1e-1;
p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%---------------------------------------------------
p.V0 = 4.5;
p.phat = 394;
p.beta = 0.3;
p.d_I = 0.1;
%p.d_V = 8.4;
p.d_V_spec = 0;
p.del_V_MPhi = 76800/200;
p.del_V_N = 2304/2.5;
%---------------------------------------------------

p.eps_L_T = 1.5*1e-5;
p.p_T_I = 0.008*2;
p.del_I_T = 238*0.5;

p = Homeostasis_calculations(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end
%----------------------------------------------------------------------

%% MODEL WITH VIRUS
tspan = [0 30];

[timeM,solM] = COVID_IMMUNE_MODEL(p,tspan);

%full_model_plotter(timeM,solM,p)
% 
%Parameters to delay IFN/get severe dynamics
p.p_F_I = 0.002;%0.0002
%p.eta_F_I = 1;%16
p.eta_F_MPhi = 1e-4*2;%1e-3
p.p_M_I = 2*0.6;

[timeS,solS] = COVID_IMMUNE_MODEL(p,tspan);

for ii = 1:length(timeM)-1
    timeM_delta(ii) = timeM(ii+1)-timeM(ii);
end
for jj = 1:length(timeS)-1
    timeS_delta(jj) = timeS(jj+1)-timeS(ii);
end

%AUCmacs_M = 1e9*sum(solM(7,1:end-1).*timeM_delta) 
%AUCmacs_S = 1e9*sum(solS(7,1:end-1).*timeS_delta)
%AUCneuts_M = 1e9*sum(solM(9,1:end-1).*timeM_delta) 
%AUCneuts_S = 1e9*sum(solS(9,1:end-1).*timeS_delta)

time_deval = linspace(tspan(1),tspan(2),1e3);
solM_deval = interp1(timeM,solM',time_deval);
solS_deval = interp1(timeS,solS',time_deval);

mild_vs_severe_plotter_single_variable2(time_deval,solM_deval',time_deval,solS_deval',p);

mild_vs_severe_with_data(time_deval,solM_deval',time_deval,solS_deval',p);

