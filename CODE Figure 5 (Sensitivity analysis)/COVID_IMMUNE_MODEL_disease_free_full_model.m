function [time, sol] = COVID_IMMUNE_MODEL_disease_free_full_model(p,tspan)
%full innate model without T cells, as of 5/18/20 

initial = [p.S0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0];
sol = ode15s(@odefun, tspan, initial);
time = linspace(tspan(1),tspan(end),1000);
solstruc = deval(sol,time);

%------------------------------------------------------------------------
function dydt = odefun(t,y,Z)

        
S = y(1);
MPhi_R = y(2);
MPhi_I = y(3);
M = y(4);
N = y(5);
T = y(6);

L_U = y(7);
L_B = y(8);
G_U = y(9);
G_B = y(10);
C_U = y(11);
C_B = y(12);
F_U = y(13);
F_B = y(14);

C_BF=C_B/(p.A_C*N);

dS = p.lam_S*(1-S/p.Smax)*S;
dMPhi_R = (1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/p.eps_V_MPhi-p.d_MPhi_R*MPhi_R;
dMPhi_I = p.p_MPhi_I_G*G_B^p.h_M_MPhi*M/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)+p.p_MPhi_I_L*L_B*M/(L_B+p.eps_L_MPhi)-(1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/p.eps_V_MPhi-p.d_MPhi_I*MPhi_I;
dM = (p.M_prod_star+(p.psi_M_max-p.M_prod_star)*G_B^p.h_M/(G_B^p.h_M+p.eps_G_M^p.h_M))*p.MR-p.p_MPhi_I_G*G_B^p.h_M_MPhi*M/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)-p.p_MPhi_I_L*L_B*M/(L_B+p.eps_L_MPhi)-p.d_M*M;
dN = (p.N_prod_star+(p.psi_N_max-p.N_prod_star)*(C_BF-p.C_BF_star)/(C_BF-p.C_BF_star+p.eps_C_N))*p.NR+p.p_N_L*L_B/(L_B+p.eps_L_N)-p.d_N*N;
%dT = p.T_prod_star-p.d_T*T;
%dT = p.T_M_prod_star+p.p_T_L*L_B*T/(L_B+p.eps_L_T)+p.p_T_F*F_B*T/(F_B+p.eps_F_T)-p.d_T*T;
dT = p.T_M_prod_star+p.p_T_F*F_B*T/(F_B+p.eps_F_T)-p.d_T*T;

dL_U = p.p_L_MPhi*MPhi_I/(MPhi_I+p.eta_L_MPhi)+p.p_L_M*M/(M+p.eta_L_M)-p.k_lin_L*L_U-p.k_B_L*((N+T+M)*p.A_L-L_B)*L_U+p.k_U_L*L_B;
dL_B = -p.k_int_L*L_B+p.k_B_L*((N+T+M)*p.A_L-L_B)*L_U-p.k_U_L*L_B;
dG_U = p.p_G_MPhi_I*MPhi_I/(MPhi_I+p.eta_G_MPhi)+p.p_G_M*M/(M+p.eta_G_M)-p.k_lin_G*G_U-p.k_B_G*(M*p.A_G-G_B)*G_U+p.k_U_G*G_B;
dG_B = -p.k_int_G*G_B+p.k_B_G*(M*p.A_G-G_B)*G_U-p.k_U_G*G_B;
dC_U = p.p_C_M*M/(M+p.eta_C_M)-p.k_lin_C*C_U-p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C+p.k_U_C*C_B;
dC_B = -p.k_int_C*C_B+p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C-p.k_U_C*C_B;
dF_U = p.p_F_MPhi*MPhi_I/(MPhi_I+p.eta_F_MPhi)+p.p_F_M*M/(M+p.eta_F_M)-p.k_lin_F*F_U-p.k_B_F*(T*p.A_F-F_B)*F_U+p.k_U_F*F_B;
dF_B = -p.k_int_F*F_B+p.k_B_F*(T*p.A_F-F_B)*F_U-p.k_U_F*F_B;

dydt = [dS;dMPhi_R;dMPhi_I;dM;dN;dT;dL_U;dL_B;dG_U;dG_B;dC_U;dC_B;dF_U;dF_B];

end
%------------------------------------------------------------------------=
%----------------------------------------------------------------------
end