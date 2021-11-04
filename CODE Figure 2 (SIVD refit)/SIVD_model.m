function [sol] = SIVD_model(p,time)

lags = [p.tau_I];
sol = dde23(@ddefun, lags, @history, time);
time = linspace(time(1),time(end),1000);
solstruc = deval(sol,time);

%------------------------------------------------------------------------
function dydt = ddefun(t,y,Z)
    ylag = Z(:,1);

    V = y(1);
    S = y(2);
    I = y(3);
    D = y(4);

    dV = p.p*I-p.d_V*V;
    dS = p.lam_S*(1-(S+I+D)/p.Smax)*S-p.beta*S*V;
    dI = p.beta*ylag(1)*ylag(2)-p.d_I*I;
    dD = p.d_I*I-p.d_D*D;
    
    dydt = [dV;dS;dI;dD];

end
%------------------------------------------------------------------------=
function s = history(t)
  s = [p.V0;p.S0;p.I0;p.D0];
end
%-------------------------------------------------------------------------
end