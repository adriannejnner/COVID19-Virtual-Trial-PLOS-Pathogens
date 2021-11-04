function [t_av, sol_av, sol_S, sol_I, sol_D, param_fit] = fit_SIVD_parameters(p,data)

%data censoring
    locs_censor = find(data.data_s5<2.5);
    data.data_s5(locs_censor) = 0;
    locs_censor = find(data.data_s6<2.5);
    data.data_s6(locs_censor) = 0;
    locs_censor = find(data.data_s18<2.5);
    data.data_s18(locs_censor) = 0;
    locs_censor = find(data.data_g1<2.5);
    data.data_g1(locs_censor) = 0;
    locs_censor = find(data.data_g2<2.5);
    data.data_g2(locs_censor) = 0;
    locs_censor = find(data.data_g5<2.5);
    data.data_g5(locs_censor) = 0;
    locs_censor = find(data.data_g6<2.5);
    data.data_g6(locs_censor) = 0;
    locs_censor = find(data.data_g7<2.5);
    data.data_g7(locs_censor) = 0;

    param_guess = [p.beta p.d_V p.d_I p.V0];% p.S0]; %[1 0.01 1e-6 0 0], [1e9 40 5 2 20]

    options = optimoptions(@lsqnonlin,'Algorithm', 'trust-region-reflective','MaxFunEval',5000);
    [param_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess,[0.01 5 0.0001 0.001],[1 20 1 12], options);   % Invoke optimizer
    
    ci = nlparci(param_fit,residual,'jacobian',jacobian)
    
    
     alpha = 0.05; %significance level
 df = length(residual) - numel(param_fit); %degrees of freedom
 crit = tinv(1-alpha/2,df);       %critical value
 covm = inv(jacobian'*jacobian) * var(residual); %covariance matrix
 [~, warnId] = lastwarn; %detect if inv threw 'nearly singular' warning.
 covmIdx = sub2ind(size(covm),1:size(covm,1),1:size(covm,2));  %indices of the diag of covm
 CI = nan(numel(param_fit),2);
 if ~strcmp(warnId, 'MATLAB:nearlySingularMatrix')
     CI(:,1) = param_fit - crit * sqrt(covm(covmIdx));
     CI(:,2) = param_fit + crit * sqrt(covm(covmIdx));
 end
 
    p.beta = param_fit(1);
   p.d_V = param_fit(2);
    p.d_I = param_fit(3);
    p.V0 = param_fit(4);
    
     alpha = 0.05; %significance level
 df = length(residual) - numel(estParams); %degrees of freedom
 crit = tinv(1-alpha/2,df);       %critical value
 covm = inv(jacobian'*jacobian) * var(residual); %covariance matrix
 [~, warnId] = lastwarn; %detect if inv threw 'nearly singular' warning.
 covmIdx = sub2ind(size(covm),1:size(covm,1),1:size(covm,2));  %indices of the diag of covm
 CI = nan(numel(estParams),2);
 if ~strcmp(warnId, 'MATLAB:nearlySingularMatrix')
     CI(:,1) = estParams - crit * sqrt(covm(covmIdx));
     CI(:,2) = estParams + crit * sqrt(covm(covmIdx));
 end
    
    [t_av, sol_av, sol_S, sol_I, sol_D] = simulation_SIVD(p,[0 p.max_time]);
%------------------------------------------------------------------------
function val = residualsfunction(param)

    p.beta = param(1);
    p.d_V = param(2);
    p.d_I = param(3);
    p.V0 = param(4);
    
    %simulate model
    [t, sol_int, sol_S, sol_I, sol_D, sol] = simulation_SIVD(p,[0 max(data.time_s5+p.lag_s5)]);
        
    %calculate residual to data
        vals5 =  [data.data_s5 -  interp1(sol.x,sol.y(1,:),data.time_s5+p.lag_s5)];
        vals6 =  [data.data_s6 -  interp1(sol.x,sol.y(1,:),data.time_s6+p.lag_s6)];
        vals18 = [data.data_s18 - interp1(sol.x,sol.y(1,:),data.time_s18+p.lag_s18)];
        valg1 =  [data.data_g1 -  interp1(sol.x,sol.y(1,:),data.time_g1+p.lag_g1)];
        valg2 =  [data.data_g2 -  interp1(sol.x,sol.y(1,:),data.time_g2+p.lag_g2)];
        valg5 =  [data.data_g5 -  interp1(sol.x,sol.y(1,:),data.time_g5+p.lag_g5)];
        valg6 =  [data.data_g6 -  interp1(sol.x,sol.y(1,:),data.time_g6+p.lag_g6)];
        valg7 =  [data.data_g7 -  interp1(sol.x,sol.y(1,:),data.time_g7+p.lag_g7)];
    %add condition for wanting S rebound

    val = [vals5,vals6,vals18,valg1,valg2,valg5,valg6,valg7];%,(sol.y(2,end-5:end)-repmat(0.16,1,6))*100];
   
    sum(abs(val))
        
end
%------------------------------------------------------------------------
end
