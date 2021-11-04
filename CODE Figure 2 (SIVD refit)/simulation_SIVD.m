function [t_av, sol_av, sol_S, sol_I, sol_D,sol] = simulation_SIVD(p,tspan)


    lags = [p.tau_I];
    sol = dde23(@ddefun, lags, @history, [tspan(1),tspan(end)]);

    t_av = linspace(0,tspan(end),100);%sol.x;
    
    sol_av = deval(sol,t_av,1);
    sol_S = deval(sol,t_av,2);
    sol_I = deval(sol,t_av,3);
    sol_D = deval(sol,t_av,4);
%------------------------------------------------------------------------
function dydt = ddefun(t,y,Z)
    ylag = Z(:,1);

    V = y(1);
    S = y(2);
    I = y(3);
    D = y(4);
    
    dV = p.p*I-p.d_V*V;
    dS = p.r*S*(1-(S+I+D)/p.Smax)-p.beta*S*V;
    dI = p.beta*ylag(1)*ylag(2)-p.d_I*I;
    dD = p.d_I*I-p.d_D*D;
    
     
    dydt = [dV;dS;dI;dD];

end   
%------------------------------------------------------------------------=
function s = history(t)
  s = [p.V0, p.S0, p.I0, p.D0];
end
%-------------------------------------------------------------------------
end