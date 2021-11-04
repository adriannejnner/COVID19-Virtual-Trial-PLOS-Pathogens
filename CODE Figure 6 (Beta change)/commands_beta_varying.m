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

% simulating different beta values

% model parameters initialised

tspan = [0 30];

% simulate model for original 0%
[time0,sol0] = COVID_IMMUNE_MODEL(p,tspan);
disp('0%')

% simulate model for 10% increase
p.beta =0.3*1.1;
[time10,sol10] = COVID_IMMUNE_MODEL(p,tspan);
disp('10%')

% simulate model for 20% increase
p.beta = 0.3*1.2;
[time20,sol20] = COVID_IMMUNE_MODEL(p,tspan);
disp('20%')

% simulate model for 30% increase
p.beta = 0.3*1.3;
[time30,sol30] = COVID_IMMUNE_MODEL(p,tspan);
disp('30%')

% simulate model for 40% increase
p.beta = 0.3*1.4;
[time40,sol40] = COVID_IMMUNE_MODEL(p,tspan);
disp('40%')

% simulate model for 50% increase
p.beta = 0.3*1.5;
[time50,sol50] = COVID_IMMUNE_MODEL(p,tspan);
disp('50%')

beta_plotted(time0,sol0,time10,sol10,time20,sol20,time30,sol30,time40,sol40,time50,sol50,p)

