function [sol_s5,sol_s6,sol_s18,sol_g1,sol_g2,sol_g5,sol_g6,sol_g7,param_fit,residual] = fitting_human_data_SIVD_lag_simul(time_mat,viral_load_mat,p)

time_s5 = time_mat{1};
time_s6 = time_mat{2};
time_s18 = time_mat{3};
time_g1 = time_mat{4};
time_g2 = time_mat{5};
time_g5 = time_mat{6};
time_g6 = time_mat{7};
time_g7 = time_mat{8};

viral_load_s5 = viral_load_mat{1};
viral_load_s6 = viral_load_mat{2};
viral_load_s18 = viral_load_mat{3};
viral_load_g1 = viral_load_mat{4};
viral_load_g2 = viral_load_mat{5};
viral_load_g5 = viral_load_mat{6};
viral_load_g6 = viral_load_mat{7};
viral_load_g7 = viral_load_mat{8};

param_guess = [p.d_V p.d_I p.beta,...
    p.lag,p.lag,p.lag,p.lag,p.lag,p.lag,p.lag,p.lag,...
    p.p, p.V0]; %p.p, p.V0,p.p, p.V0,p.p, p.V0,p.p, p.V0,...

lb = [0 0 0,...
    0 0 0 0 0 0 0 0,...
    1 0.00001];
   % 1 0.00001 1 0.00001 1 0.00001 1 0.00001,...
ub = [50 5 5,...
    9 9 9 9 9 9 9 9,...
    1e3 12 ];
   % 1e3 12 1e3 12 1e3 12,...

options = optimoptions(@lsqnonlin,'Algorithm', 'trust-region-reflective','MaxFunEval',5000);
residualsfunction2 = @(param)residualsfunction(param);
[param_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(residualsfunction2, param_guess,lb,ub,options);   % Invoke optimizer


    p.d_V = param_fit(1);
    p.d_I = param_fit(2);
    p.beta = param_fit(3);
    
    p.lag_s5 = param_fit(4);
    p.lag_s6 = param_fit(5);
    p.lag_s18 = param_fit(6);
    p.lag_g1 = param_fit(7);
    p.lag_g2 = param_fit(8);
    p.lag_g5 = param_fit(9);
    p.lag_g6 = param_fit(10);
    p.lag_g7 = param_fit(11);
    
    p.p   = param_fit(12);
    
    p.V0  = param_fit(13);
    [sol_s5] = SIVD_model(p,[0 time_s5(end)+p.lag_s5]);
   % p.p   = param_fit(14);
   % p.V0  = param_fit(14);
    [sol_s6] = SIVD_model(p,[0 time_s6(end)+p.lag_s6]);
   % p.p   = param_fit(16);
   % p.V0  = param_fit(15);
    [sol_s18] = SIVD_model(p,[0 time_s18(end)+p.lag_s18]);
   % p.p   = param_fit(18);
   % p.V0  = param_fit(16);
    [sol_g1] = SIVD_model(p,[0 time_g1(end)+p.lag_g1]);
   % p.p   = param_fit(20);
   % p.V0  = param_fit(17);
    [sol_g2] = SIVD_model(p,[0 time_g2(end)+p.lag_g2]);
   % p.p   = param_fit(22);
   % p.V0  = param_fit(18);
    [sol_g5] = SIVD_model(p,[0 time_g5(end)+p.lag_g5]);
   % p.p   = param_fit(24);
   % p.V0  = param_fit(19);
    [sol_g6] = SIVD_model(p,[0 time_g6(end)+p.lag_g6]);
   % p.p   = param_fit(26);
   % p.V0  = param_fit(20);
    [sol_g7] = SIVD_model(p,[0 time_g7(end)+p.lag_g7]);

%------------------------------------------------------------------------
function val = residualsfunction(param)

    p.d_V = param(1);
    p.d_I = param(2);
    p.beta = param(3);
    
    p.lag_s5 = param(4);
    p.lag_s6 = param(5);
    p.lag_s18 = param(6);
    p.lag_g1 = param(7);
    p.lag_g2 = param(8);
    p.lag_g5 = param(9);
    p.lag_g6 = param(10);
    p.lag_g7 = param(11);
    
    p.p   = param(12);
    p.V0  = param(13);
    [sol] = SIVD_model(p,[0 time_s5(end)+p.lag_s5]);
     val_s5 = [viral_load_s5 - deval(sol,time_s5+p.lag_s5,1)];
     
  %  p.p   = param(14);
  %  p.V0  = param(14);
    [sol] = SIVD_model(p,[0 time_s6(end)+p.lag_s6]);
     val_s6 = [viral_load_s6 - deval(sol,time_s6+p.lag_s6,1)];
     
  %  p.p   = param(16);
  %  p.V0  = param(15);
    [sol] = SIVD_model(p,[0 time_s18(end)+p.lag_s18]);
     val_s18 = [viral_load_s18 - deval(sol,time_s18+p.lag_s18,1)];
     
  %  p.p   = param(18);
  %  p.V0  = param(16);
    [sol] = SIVD_model(p,[0 time_g1(end)+p.lag_g1]);
     val_g1 = [viral_load_g1 - deval(sol,time_g1+p.lag_g1,1)];
    
  %  p.p   = param(20);
  %  p.V0  = param(17);
    [sol] = SIVD_model(p,[0 time_g2(end)+p.lag_g2]);
     val_g2 = [viral_load_g2 - deval(sol,time_g2+p.lag_g2,1)];
    
  %  p.p   = param(22);
  %  p.V0  = param(18);     
    [sol] = SIVD_model(p,[0 time_g5(end)+p.lag_g5]);
     val_g5 = [viral_load_g5 - deval(sol,time_g5+p.lag_g5,1)];
     
   % p.p   = param(24);
   % p.V0  = param(19);
    [sol] = SIVD_model(p,[0 time_g6(end)+p.lag_g6]);
     val_g6 = [viral_load_g6 - deval(sol,time_g6+p.lag_g6,1)];
     
   % p.p   = param(26);
   % p.V0  = param(20);
    [sol] = SIVD_model(p,[0 time_g7(end)+p.lag_g7]);
     val_g7 = [viral_load_g7 - deval(sol,time_g7+p.lag_g7,1)];

     val = [val_s5,val_s6,val_s18,val_g1,val_g2,val_g5,val_g6,val_g7];
     sum(abs(val))
  
end
%------------------------------------------------------------------------
end