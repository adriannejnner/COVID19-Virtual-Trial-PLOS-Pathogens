function p = Homeostasis_calculations(p)

% << RESIDENT MACROPHAGES >>
% assume that at homeostasis there is no death of resident macrophages, or
% if there is, it is negligable
p.d_MPhi_R = 0;

% <<MONOCYTES>>
p.M_star = p.M0;  % homeostasis for monocytes
p.MPhi_R_star = p.MPhi_R_0;  % homeostasis for monocytes

% <<NEUTROPHILS>>
p.N_star = p.N0;           % homeostatic concentration of neutrophils

% << T cells >>
p.T_star = p.T0;

% << IFN >>
p.F_U_star = p.F_U_0;
p.F_B_star = p.k_B_F*p.T_star.*p.A_F*p.F_U_star/(p.k_int_F+p.k_B_F*p.F_U_star+p.k_U_F);
p.F_B_0 = p.F_B_star;

% << GM-CSF>>
p.G_U_star = p.G_U_0;    % homeostatic unbound concentration
p.G_B_star = (p.k_B_G*p.M_star*p.A_G*p.G_U_star)/(p.k_int_G+p.k_U_G+p.k_B_G*p.G_U_star);  % homeostatic bound concentration
p.G_B_0 = p.G_B_star;       % Initial G-CSF (bound) in ng/mL

% << G-CSF>>
p.C_U_star = p.C_U_0;    % homeostatic unbound concentration
p.C_B_star=(p.k_B_C*p.C_U_star^p.stoch_C*p.A_C*(p.N_star))/(p.k_int_C+p.k_U_C+p.C_U_star^p.stoch_C*p.k_B_C);
p.C_B_0 = p.C_B_star;       % Initial G-CSF (bound) in ng/mL

%p.C_prod_star=p.k_lin_C*p.C_U_star+p.k_B_C*p.C_U_star^p.stoch_C*p.A_C*(p.N_star)-p.k_B_C*p.C_U_star^p.stoch_C*p.C_B_star-p.k_U_C*p.C_B_star;
p.C_BF_star=p.C_B_star/(p.A_C*p.N_star); % homeostatic bound fraction (unitless)

% <<IL-6>>
p.L_U_star = p.L_U_0;
p.L_B_star = p.k_B_L*(p.T_star+p.N_star+p.M_star)*p.A_L*p.L_U_star/(p.k_int_L+p.k_B_L*p.L_U_star+p.k_U_L);
p.L_B_0 = p.L_B_star;

% <<SYSTEMIC MACROPHAGES>>
p.MPhi_I_star = (p.p_MPhi_I_G*p.G_B_star^p.h_M_MPhi*p.M_star/(p.G_B_star^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)+p.p_MPhi_I_L*p.L_B_star*p.M_star/(p.L_B_star+p.eps_L_MPhi))/((1-p.MPhi_R_star/p.MPhi_max)*p.lam_MPhi/p.eps_V_MPhi+p.d_MPhi_I);
p.MPhi_I_0 = p.MPhi_I_star;

% << IL-6 >>
% estimating production of IL-6 by macrophages, assuming that in the
% homeostastic system this may be different
%p.p_L_MPhi = (p.MPhi_I_star+p.eta_L_MPhi)/p.MPhi_I_star*(-(p.p_L_M*p.M_star/(p.M_star+p.eta_L_M)-p.k_lin_L*p.L_U_star-p.k_B_L*((p.N_star+p.T_star+p.M_star)*p.A_L-p.L_B_star)*p.L_U_star+p.k_U_L*p.L_B_star));
KK = 1/p.p_L_MPhi*(-p.p_L_M*p.M_star/(p.M_star+p.eta_L_M)+p.k_lin_L*p.L_U_star+p.k_B_L*((p.N_star+p.T_star+p.M_star)*p.A_L-p.L_B_star)*p.L_U_star-p.k_U_L*p.L_B_star);
p.eta_L_MPhi = (p.MPhi_I_star-KK*p.MPhi_I_star)/KK;

% << GM-CSF >>
p.eta_G_MPhi = p.eta_L_MPhi;
p.p_G_MPhi_I = -(-p.k_lin_G*p.G_U_star-p.k_B_G*(p.M_star*p.A_G-p.G_B_star)*p.G_U_star+p.k_U_G*p.G_B_star)/(p.MPhi_I_star/(p.MPhi_I_star+p.eta_G_MPhi)+p.M_star/(p.M_star+p.eta_G_M));
p.p_G_M = p.p_G_MPhi_I; %this was assuming that p_G_MPhi_I = p_G_M
p.p_C_M = p.p_G_M/100; %unknown             % production rate of G-CSF by monocyte XXXXXXXXXXXXXX                           UNKNOWN

%estimating half-effect of Monocyte production of G-CSF
BB = p.k_lin_C*p.C_U_star+p.k_B_C*(p.N_star*p.A_C-p.C_B_star)*(p.C_U_star).^p.stoch_C-p.k_U_C*p.C_B_star;
p.eta_C_M = (p.p_C_M*p.M_star-p.M_star*BB)/BB;

% << T cells >>
p.T_prod_star = p.d_T*p.T_star;
%p.T_M_prod_star = p.d_T*p.T_star-p.p_T_L*p.L_B_star*p.T_star/(p.L_B_star+p.eps_L_T)-p.p_T_F*p.F_B_star*p.T_star/(p.F_B_star+p.eps_F_T);
p.T_M_prod_star = p.d_T*p.T_star-p.p_T_F*p.F_B_star*p.T_star/(p.F_B_star+p.eps_F_T);

% <<MONOCYTES>>
DD = 1/p.MR*(p.p_MPhi_I_G*p.G_B_star^p.h_M_MPhi*p.M_star/(p.G_B_star^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)+p.p_MPhi_I_L*p.L_B_star*p.M_star/(p.L_B_star+p.eps_L_MPhi)+p.d_M*p.M_star);
EE = p.G_B_star^p.h_M/(p.G_B_star^p.h_M+p.eps_G_M^p.h_M);
p.M_prod_star = (DD-p.psi_M_max*EE)/(1-EE);

% <<NEUTROPHILS>>
%p.N_prod_star = (p.d_N*p.N_star-p.p_N_L*p.L_B_star/(p.L_B_star+p.eps_L_N))*1/p.NR;  % homeostasis production rate for neutrophils
p.N_prod_star = p.d_N*p.N_star/(p.NR+p.L_B_star/(p.L_B_star+p.eps_L_N));
p.p_N_L = p.N_prod_star;

% << IFN >>
% estimating half effect of macrophage production of IFN assuming that
% macrophages and monocytes produce IFN at the same rate
%p.p_F_MPhi = p.p_F_M;
CC = p.p_F_M*p.M_star/(p.M_star+p.eta_F_M)-p.k_lin_F*p.F_U_star-p.k_B_F*(p.T_star*p.A_F-p.F_B_star)*p.F_U_star+p.k_U_F*p.F_B_star;
p.eta_F_MPhi = (p.p_F_MPhi*p.MPhi_I_star+CC*p.MPhi_I_star)/(-CC);


p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;

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

p.p_G_M = 1.234*1e3*1e-1;
%p.p_G_MPhi_I = 1232.4;

p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

end