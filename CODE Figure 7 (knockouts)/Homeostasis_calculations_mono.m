function p = Homeostasis_calculations(p)

% << RESIDENT MACROPHAGES >>
% assume that at homeostasis there is no death of resident macrophages, or
% if there is, it is negligable
p.d_MPhi_R = 0;

% <<MONOCYTES>>
p.M_star = 0;  % homeostasis for monocytes
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

end