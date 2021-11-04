function [time, sol] = COVID_IMMUNE_MODEL_virus_IFN_resistance(p,tspan)

opts = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
solstruc = ddesd_f5(@ddefun,@(t,y) delayP(t,y,p),@history,tspan,opts);

time = linspace(tspan(1),tspan(end),1000);
sol = deval(solstruc,time);

%------------------------------------------------------------------------
function dydt = ddefun(t,y,Z)
ylag1 = Z(:,1);
ylag2 = Z(:,2);
        
V = y(1);
S = y(2);
I = y(3);
R = y(4);
D = y(5);
F_U = y(6);
F_B = y(7);
      

if t>1
if I<1e-9
    I = 0;
end
if D<1e-9
    D = 0;
end
end

dV = p.phat*I-p.d_V*V;%/(1+F_B/p.eps_F_I)
dS = p.lam_S*(1-(S+I+D+R)/p.Smax)*S-p.beta*S*V;%/(1+F_B/p.eps_F_I);%/(1+F_B/p.eps_F_I)
dI =p.beta/(1+F_B/p.eps_F_I)*ylag1(2)*ylag1(1)-p.d_I*I;%
dR =p.lam_S*(1-(S+I+D+R)/p.Smax)*R+p.beta*S*V/(1+p.eps_F_I/F_B);%-0.001*R;
dD = p.d_I*I-p.d_D*D;

dF_U =p.psi_F_prod+p.p_F_I*I/(I+p.eta_F_I)-p.k_lin_F*F_U-p.k_B_F*((p.T_star+I)*p.A_F-F_B)*F_U+p.k_U_F*F_B;
dF_B = -p.k_int_F*F_B+p.k_B_F*((p.T_star+I)*p.A_F-F_B)*F_U-p.k_U_F*F_B;

[p.psi_F_prod p.p_F_I*I/(I+p.eta_F_I)];

if t>1
if I<1e-9
    dI = 0;
end
if D<1e-9
    dD = 0;
end
end
dydt = [dV;dS;dI;dR;dD;dF_U;dF_B];
end
%------------------------------------------------------------------------=
function s = history(t)
  s = [p.V0;p.S0;p.I0;0;p.D0;p.F_U_0;p.F_B_0];
end
function d = delayP(t,y,p)
%This function sets up the delay vectors necessary for the DDE solver.
d = [t-p.tau_I,t-p.tau_T];     
end
%-------------------------------------------------------------------------
end