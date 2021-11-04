function [time,sol,solstruc] = COVID_IMMUNE_MODEL(p,tspan)

opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
solstruc = ddesd_f5(@ddefun,@(t,y) delayP(t,y,p),@history,tspan,opts);

time = solstruc.x;%linspace(tspan(1),tspan(end),1000);
sol = solstruc.y;%deval(solstruc,time);

%------------------------------------------------------------------------
function dydt = ddefun(t,y,Z)
ylag1 = Z(:,1);
ylag2 = Z(:,2);
    
V = y(1);
S = y(2);
I = y(3);
R = y(4);
D = y(5);
MPhi_R = y(6);
MPhi_I = y(7);
M = y(8);
N = y(9);
T = y(10);

L_U = y(11);
L_B = y(12);
G_U = y(13);
G_B = y(14);
C_U = y(15);
C_B = y(16);
F_U = y(17);
F_B = y(18);

C_BF=C_B/(p.A_C*N);
        
% 
% if t>8
%     
% if I<1e-9;
%       I = 0;
%  end
% if D<1e-9;
%       D = 0;
% end
% end
% 
% if t>8
% if    V<1e-10
%     V = 0;
% end
% end

if isempty(find(y<0))==0
   disp('someone negative') 
end

MPhi_I = 0;

dV = p.phat*I-p.del_V_MPhi*MPhi_I*V-p.del_V_N*V*N-p.d_V_spec*V;
dS = p.lam_S*(1-(S+I+D+R)/p.Smax)*S-p.beta*S*V-p.rho*p.del_N/(1+(p.IC_50_N/N)^p.h_N)*S;%-p.rho*p.del_N/(1+(p.IC_50_N/M)^p.h_N)*S-p.rho*p.del_N/(1+(0.009/MPhi_I)^p.h_N)*S;
dI = p.beta/(1+F_B/p.eps_F_I)*ylag1(2)*ylag1(1)-p.d_I*I-p.del_N/(1+(p.IC_50_N/N).^p.h_N)*I-p.del_I_MPhi*MPhi_I*I-p.del_I_T*T*I;   
dR = p.lam_S*(1-(S+I+D+R)/p.Smax)*R+p.beta*S*V/(1+p.eps_F_I/F_B)-p.rho*p.del_N/(1+(p.IC_50_N*2/N)^p.h_N)*R+p.lam_S*(1-(S+I+D+R)/p.Smax)*R;
dD = p.d_I*I+p.del_N*I/(1+(p.IC_50_N/N)^p.h_N)+p.del_N*p.rho*S/(1+(p.IC_50_N/N)^p.h_N)+p.del_I_MPhi*MPhi_I*I+p.del_I_T*T*I-p.d_D*D-p.del_D_MPhi*D*(MPhi_R+MPhi_I)+p.del_MPhi_D*D*(MPhi_I)+p.del_N*p.rho/(1+(p.IC_50_N/N)^p.h_N)*R;%+p.rho*p.del_N/(1+(p.IC_50_N/M)^p.h_N)*S+p.rho*p.del_N/(1+(0.009/MPhi_I)^p.h_N)*S;%close+p.d_M*M+p.d_N*N+p.d_T*T;

dMPhi_R = 0;%-p.a_I_MPhi*MPhi_R*(I+D)+(1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/(V+p.eps_V_MPhi);
dMPhi_I = 0;%p.a_I_MPhi*MPhi_R*(I+D)+p.p_MPhi_I_G*G_B^p.h_M_MPhi/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)*M+p.p_MPhi_I_L*L_B/(L_B+p.eps_L_MPhi)*M-p.d_MPhi_I*MPhi_I-p.del_MPhi_D*MPhi_I*D-(1-MPhi_R/p.MPhi_max)*p.lam_MPhi*MPhi_I/(V+p.eps_V_MPhi);%-p.del_I_T*T*MPhi_I*100;%-p.del_MPhi_D*MPhi_I*T;    
dM = (p.M_prod_star+(p.psi_M_max-p.M_prod_star)*G_B^p.h_M/(G_B^p.h_M+p.eps_G_M^p.h_M))*p.MR+p.p_M_I*I*M/(I+p.eps_I_M)-p.p_MPhi_I_G*G_B^p.h_M_MPhi*M/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)-p.p_MPhi_I_L*L_B*M/(L_B+p.eps_L_MPhi)-p.d_M*M;
dN = (p.N_prod_star+(p.psi_N_max-p.N_prod_star)*(C_BF-p.C_BF_star)/(C_BF-p.C_BF_star+p.eps_C_N))*p.NR+p.p_N_L*L_B/(L_B+p.eps_L_N)-p.d_N*N;
dT = p.p_T_I*ylag2(3)/(1+L_B/p.eps_L_T)+p.p_T_F*F_B/(F_B+p.eps_F_T)*T-p.d_T*T;

%[t p.a_I_MPhi*MPhi_R*(I+D) p.p_MPhi_I_G*G_B^p.h_M_MPhi/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)*M p.p_MPhi_I_L*L_B/(L_B+p.eps_L_MPhi)*M]
%[t p.p_MPhi_I_G*G_B^p.h_M_MPhi/(G_B^p.h_M_MPhi+(p.eps_G_MPhi/2)^p.h_M_MPhi)*M p.p_MPhi_I_G*G_B^p.h_M_MPhi/(G_B^p.h_M_MPhi+p.eps_G_MPhi^p.h_M_MPhi)*M  p.p_MPhi_I_L*L_B/(L_B+p.eps_L_MPhi)*M]

dL_U = p.p_L_I*I/(I+p.eta_L_I)+p.p_L_MPhi*MPhi_I/(MPhi_I+p.eta_L_MPhi)+p.p_L_M*M/(M+p.eta_L_M)-p.k_lin_L*L_U-p.k_B_L*((N+T+M)*p.A_L-L_B)*L_U+p.k_U_L*L_B;
dL_B = -p.k_int_L*L_B+p.k_B_L*((N+T+M)*p.A_L-L_B)*L_U-p.k_U_L*L_B;
dG_U = p.p_G_MPhi_I*MPhi_I/(MPhi_I+p.eta_G_MPhi)+p.p_G_M*M/(M+p.eta_G_M)-p.k_lin_G*G_U-p.k_B_G*(M*p.A_G-G_B)*G_U+p.k_U_G*G_B;
dG_B = -p.k_int_G*G_B+p.k_B_G*(M*p.A_G-G_B)*G_U-p.k_U_G*G_B;
dC_U = p.p_C_M*M/(M+p.eta_C_M)-p.k_lin_C*C_U-p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C+p.k_U_C*C_B;
dC_B = -p.k_int_C*C_B+p.k_B_C*(N*p.A_C-C_B)*C_U^p.stoch_C-p.k_U_C*C_B;
dF_U = p.p_F_I*I/(I+p.eta_F_I)+p.p_F_MPhi*MPhi_I/(MPhi_I+p.eta_F_MPhi)+p.p_F_M*M/(M+p.eta_F_M)-p.k_lin_F*F_U-p.k_B_F*((T+I)*p.A_F-F_B)*F_U+p.k_U_F*F_B;
dF_B = -p.k_int_F*F_B+p.k_B_F*((T+I)*p.A_F-F_B)*F_U-p.k_U_F*F_B;
% 
% if t>8
%     
% if I<1e-9;
%       dI = 0;
%  end
% if D<1e-9;
%       dD = 0;
% end
% end
dydt = [dV;dS;dI;dR;dD;dMPhi_R;dMPhi_I;dM;dN;dT;dL_U;dL_B;dG_U;dG_B;dC_U;dC_B;dF_U;dF_B];

end
%------------------------------------------------------------------------=
function s = history(t)
  s = [p.V0;p.S0;p.I0;p.R0;p.D0;p.MPhi_R_0;p.MPhi_I_0;p.M0;p.N0;p.T0;p.L_U_0;p.L_B_0;p.G_U_0;p.G_B_0;p.C_U_0;p.C_B_0;p.F_U_0;p.F_B_0];
end
function d = delayP(t,y,p)
%This function sets up the delay vectors necessary for the DDE solver.
d = [t-p.tau_I,t-p.tau_T];     
end
%-------------------------------------------------------------------------
end