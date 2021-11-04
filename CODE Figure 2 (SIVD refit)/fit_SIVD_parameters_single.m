function [t_av, sol_av, sol_S, sol_I, sol_D, param_fit] = fit_SIVD_parameters(p,data,time)

%data censoring
    locs_censor = find(data<2.5);
    data(locs_censor) = 0;

    param_guess = [p.beta p.d_V p.d_I p.V0];% p.S0]; %[1 0.01 1e-6 0 0], [1e9 40 5 2 20]

    options = optimoptions(@lsqnonlin,'Algorithm', 'trust-region-reflective','MaxFunEval',5000);
    [param_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess,[0.01 5 0.0001 0.001],[1 20 1 12], options);   % Invoke optimizer

    p.beta = param_fit(1);
   p.d_V = param_fit(2);
    p.d_I = param_fit(3);
    p.V0 = param_fit(4);
    
    [t_av, sol_av, sol_S, sol_I, sol_D] = simulation_SIVD(p,[0 p.max_time]);
%------------------------------------------------------------------------
function val = residualsfunction(param)

    p.beta = param(1);
    p.d_V = param(2);
    p.d_I = param(3);
    p.V0 = param(4);
    
    %simulate model
    [t, sol_int, sol_S, sol_I, sol_D, sol] = simulation_SIVD(p,[0 max(time)]);
        
    %calculate residual to data
        val =  [data -  interp1(sol.x,sol.y(1,:),time)];
     %add condition for wanting S rebound

    sum(abs(val))
        
end
%------------------------------------------------------------------------
end
