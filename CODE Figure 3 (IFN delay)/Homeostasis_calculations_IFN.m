function p = Homeostasis_calculations_IFN(p)

p.A_F = p.stoch*p.MM_F/p.avo*(p.R_F_I+p.R_F_T)*(1/5000)*10^9*1e12;% IFN

% << RESIDENT MACROPHAGES >>
% assume that at homeostasis there is no death of resident macrophages, or
% if there is, it is negligable
p.d_MPhi_R = 0;

% <<MONOCYTES>>
p.M_star =0;  % homeostasis for monocytes
p.MPhi_R_star = 0;  % homeostasis for monocytes

% <<NEUTROPHILS>>
p.N_star = 0;           % homeostatic concentration of neutrophils

% << T cells >>
p.T_star = p.T0;

% << IFN >>
p.F_U_star = p.F_U_0;
p.F_B_star = p.k_B_F*p.T_star.*p.A_F*p.F_U_star/(p.k_int_F+p.k_B_F*p.F_U_star+p.k_U_F);
p.F_B_0 = p.F_B_star;

p.psi_F_prod = p.k_lin_F*p.F_U_star+p.k_B_F*(p.T_star*p.A_F-p.F_B_star)*p.F_U_star-p.k_U_F*p.F_B_star;

end