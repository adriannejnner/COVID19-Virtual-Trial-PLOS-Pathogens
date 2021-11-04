%simulating full model

format long

p = load_parameters();

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;
p.p_F_M = 3.5600*1e4*1e-4; 
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;
p.eps_G_MPhi = 2664.5/1e5;
p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;
p.p_T_I = 0.01*0.8;
p.del_I_T = 238;
p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; 
p.eps_I_M = 0.11;
p.del_V_N = 768*3;
p.del_V_MPhi = 768*100;
p.p_G_M = 1.234*1e3*1e-1;
p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%---------------------------------------------------
p.V0 = 4.5;
p.phat = 394;
p.beta = 0.3;
p.d_I = 0.1;
%p.d_V = 8.4;
p.d_V_spec = 0;
p.del_V_MPhi = 76800/200;
p.del_V_N = 2304/2.5;
%---------------------------------------------------

p.eps_L_T = 1.5*1e-5;
p.p_T_I = 0.008*2;
p.del_I_T = 238*0.5;

p = Homeostasis_calculations(p);

%-----------------------------------------------------------------------
estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
end
%----------------------------------------------------------------------

tol = 5;
tspan = [0 30];
[timeORIG,solORIG] = COVID_IMMUNE_MODEL_FM(p,tspan);% solve model

factor = 0.2;% 0.1 - 10%

tunder =  find([solORIG(2,:)+solORIG(4,:)]<=0.16*0.3);
tabove =  find(solORIG(2,tunder:end)+solORIG(4,tunder:end)>=0.16*0.3);
if isempty(tabove)==1 & isempty(tunder)==1
   p.sol_orig_under = 0;
else
    if isempty(tabove)==1
        tabove = 30;
    end
    p.sol_orig_under = timeORIG(tabove(1)+tunder(1))-timeORIG(tunder(1));
end


for i = 1:length(solORIG(17,:))-1
   tgrid_vec(i) = timeORIG(i+1)-timeORIG(i); 
end
IFNexposure_ORIG = sum(solORIG(17,1:end-1).*tgrid_vec);
figure


peak_loc = find(solORIG(17,:)==max(solORIG(17,:)));
p.sol_orig_peak = timeORIG(peak_loc); %time of IFN peak
tol = 5;


% PHAT
p.phat = p.phat*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.phat = p.phat/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.phat = p.phat/(1-factor);%revert to old parameter
mega_metric_matrix(1:2,:) = [metric_vector_10plus;metric_vector_10minus];

% epsilon_F_I
p.eps_F_I = p.eps_F_I*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.eps_F_I = p.eps_F_I/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.eps_F_I = p.eps_F_I/(1-factor);%revert to old parameter
mega_metric_matrix(3:4,:) = [metric_vector_10plus;metric_vector_10minus];

% d_V_spec
p.d_V_spec = p.d_V_spec*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.d_V_spec = p.d_V_spec/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.d_V_spec = p.d_V_spec/(1-factor);%revert to old parameter
mega_metric_matrix(5:6,:) = [metric_vector_10plus;metric_vector_10minus];

% del_V_MPhi
p.del_V_MPhi = p.del_V_MPhi*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_V_MPhi = p.del_V_MPhi/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_V_MPhi = p.del_V_MPhi/(1-factor);%revert to old parameter
mega_metric_matrix(7:8,:) = [metric_vector_10plus;metric_vector_10minus];

% del_V_N
p.del_V_N = p.del_V_N*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_V_N = p.del_V_N/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_V_N = p.del_V_N/(1-factor);%revert to old parameter
mega_metric_matrix(9:10,:) = [metric_vector_10plus;metric_vector_10minus];

% beta
p.beta = p.beta*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.beta = p.beta/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.beta = p.beta/(1-factor);%revert to old parameter
mega_metric_matrix(11:12,:) = [metric_vector_10plus;metric_vector_10minus];

% lam_S
p.lam_S = p.lam_S*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.lam_S = p.lam_S/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.lam_S = p.lam_S/(1-factor);%revert to old parameter
mega_metric_matrix(13:14,:) = [metric_vector_10plus;metric_vector_10minus];

% delta_N
p.del_N = p.del_N*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_N = p.del_N/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_N = p.del_N/(1-factor);%revert to old parameter
mega_metric_matrix(15:16,:) = [metric_vector_10plus;metric_vector_10minus];

% delta_N
p.IC_50_N = p.IC_50_N*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.IC_50_N = p.IC_50_N/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.IC_50_N = p.IC_50_N/(1-factor);%revert to old parameter
mega_metric_matrix(17:18,:) = [metric_vector_10plus;metric_vector_10minus];

% d_I
p.d_I = p.d_I*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.d_I = p.d_I/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.d_I = p.d_I/(1-factor);%revert to old parameter
mega_metric_matrix(19:20,:) = [metric_vector_10plus;metric_vector_10minus];

% deltaIMPhi
p.del_I_MPhi = p.del_I_MPhi*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_I_MPhi = p.del_I_MPhi/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_I_MPhi = p.del_I_MPhi/(1-factor);%revert to old parameter
mega_metric_matrix(21:22,:) = [metric_vector_10plus;metric_vector_10minus];

% deltaIT
p.del_I_T = p.del_I_T*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_I_T = p.del_I_T/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_I_T = p.del_I_T/(1-factor);%revert to old parameter
mega_metric_matrix(23:24,:) = [metric_vector_10plus;metric_vector_10minus];

% delta_MPhi_D
p.del_MPhi_D = p.del_MPhi_D*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_MPhi_D = p.del_MPhi_D/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_MPhi_D = p.del_MPhi_D/(1-factor);%revert to old parameter
mega_metric_matrix(25:26,:) = [metric_vector_10plus;metric_vector_10minus];

% delta_MPhi_D
p.del_D_MPhi = p.del_D_MPhi*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_D_MPhi = p.del_D_MPhi/(1+factor)*(1-factor); %decrease by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.del_D_MPhi = p.del_D_MPhi/(1-factor);%revert to old parameter
mega_metric_matrix(27:28,:) = [metric_vector_10plus;metric_vector_10minus];

% a_I_MPhi
p.a_I_MPhi = p.a_I_MPhi*(1+factor); % increase by 10%
metric_vector_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.a_I_MPhi = p.a_I_MPhi/(1+factor)*(1-factor); % increase by 10%
metric_vector_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);% calculate all the metrics
p.a_I_MPhi = p.a_I_MPhi/(1-factor);%revert to old parameter
mega_metric_matrix(29:30,:) = [metric_vector_10plus;metric_vector_10minus];

% eps_V_MPhi
p.eps_V_MPhi = p.eps_V_MPhi*(1+factor); % increase by 10%
metric_vectorEpsVMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_V_MPhi = p.eps_V_MPhi/(1+factor)*(1-factor);
metric_vectorEpsVMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_V_MPhi = p.eps_V_MPhi/(1-factor);%revert to old parameter
mega_metric_matrix(31:32,:) = [metric_vectorEpsVMPhi_10plus;metric_vectorEpsVMPhi_10minus];

% eps_G_MPhi
p.eps_G_MPhi = p.eps_G_MPhi*(1+factor); % increase by 10%
metric_vectorEpsGMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_G_MPhi = p.eps_G_MPhi/(1+factor)*(1-factor);
metric_vectorEpsGMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_G_MPhi = p.eps_G_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix(33:34,:) = [metric_vectorEpsGMPhi_10plus;metric_vectorEpsGMPhi_10minus];

% P_g_MPhi
p.p_MPhi_I_G = p.p_MPhi_I_G*(1+factor); % increase by 10%
metric_vectorPGMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_MPhi_I_G = p.p_MPhi_I_G/(1+factor)*(1-factor);
metric_vectorPGMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_MPhi_I_G = p.p_MPhi_I_G/(1-factor);%revert to old parameter

mega_metric_matrix(35:36,:) = [metric_vectorPGMPhi_10plus;metric_vectorPGMPhi_10minus];

% eps_L_MPhi
p.eps_L_MPhi = p.eps_L_MPhi*(1+factor); % increase by 10%
metric_vectorEpsLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_MPhi = p.eps_L_MPhi/(1+factor)*(1-factor);
metric_vectorEpsLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_MPhi = p.eps_L_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix(37:38,:) = [metric_vectorEpsLMPhi_10plus;metric_vectorEpsLMPhi_10minus];

% P_L_MPhi
p.p_MPhi_I_L = p.p_MPhi_I_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_MPhi_I_L = p.p_MPhi_I_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_MPhi_I_L = p.p_MPhi_I_L/(1-factor);%revert to old parameter

mega_metric_matrix(39:40,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('39')
% eps_G_M
p.eps_G_M = p.eps_G_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_G_M = p.eps_G_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_G_M = p.eps_G_M/(1-factor);%revert to old parameter

mega_metric_matrix(41:42,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('41')

%
% psi_M_Max
p.psi_M_max = p.psi_M_max*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.psi_M_max = p.psi_M_max/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.psi_M_max = p.psi_M_max/(1-factor);%revert to old parameter

mega_metric_matrix(43:44,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('43')

% p_M_I
p.p_M_I = p.p_M_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_M_I = p.p_M_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_M_I = p.p_M_I/(1-factor);%revert to old parameter
mega_metric_matrix(45:46,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('45')

% eps_I_M
p.eps_I_M = p.eps_I_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_I_M = p.eps_I_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_I_M = p.eps_I_M/(1-factor);%revert to old parameter

mega_metric_matrix(47:48,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('47')
% psi_N_max
p.psi_N_max = p.psi_N_max*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.psi_N_max = p.psi_N_max/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.psi_N_max = p.psi_N_max/(1-factor);%revert to old parameter

mega_metric_matrix(49:50,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('49')
% eps_C_N
p.eps_C_N = p.eps_C_N*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_C_N = p.eps_C_N/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_C_N = p.eps_C_N/(1-factor);%revert to old parameter

mega_metric_matrix(51:52,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('51')
% p_N_L
p.p_N_L = p.p_N_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_N_L = p.p_N_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_N_L = p.p_N_L/(1-factor);%revert to old parameter

mega_metric_matrix(53:54,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('53')
% eps_L_N
p.eps_L_N = p.eps_L_N*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_N = p.eps_L_N/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_N = p.eps_L_N/(1-factor);%revert to old parameter

mega_metric_matrix(55:56,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('55')

% eps_L_T
p.eps_L_T = p.eps_L_T*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_T = p.eps_L_T/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_L_T = p.eps_L_T/(1-factor);%revert to old parameter

mega_metric_matrix(57:58,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('57')
% p_T_I
p.p_T_I = p.p_T_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_T_I = p.p_T_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_T_I = p.p_T_I/(1-factor);%revert to old parameter

mega_metric_matrix(59:60,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('59')
%
%eps_F_T
p.eps_F_T = p.eps_F_T*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_F_T = p.eps_F_T/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eps_F_T = p.eps_F_T/(1-factor);%revert to old parameter

mega_metric_matrix(61:62,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('61')

% p_T_F
p.p_T_F = p.p_T_F*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_T_F = p.p_T_F/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_T_F = p.p_T_F/(1-factor);%revert to old parameter

mega_metric_matrix(63:64,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('75')

% V0
p.V0 = p.V0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.V0 = p.V0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.V0 = p.V0/(1-factor);%revert to old parameter

mega_metric_matrix(65:66,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('63')
% S0
p.S0 = p.S0*(1+factor); % increase by 10%
p.Smax = p.S0;
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.S0 = p.S0/(1+factor)*(1-factor);
p.Smax = p.S0;
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.S0 = p.S0/(1-factor);%revert to old parameter
p.Smax = p.S0;

mega_metric_matrix(67:68,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('65')

% MPhiR0
p.MPhi_R_0 = p.MPhi_R_0*(1+factor); % increase by 10%
p.MPhi_max = p.MPhi_R_0;
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.MPhi_R_0 = p.MPhi_R_0/(1+factor)*(1-factor);
p.MPhi_max = p.MPhi_R_0;
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.MPhi_R_0 = p.MPhi_R_0/(1-factor);%revert to old parameter
p.MPhi_max = p.MPhi_R_0;

mega_metric_matrix(69:70,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('67')
% M0
p.M0 = p.M0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.M0 = p.M0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.M0 = p.M0/(1-factor);%revert to old parameter

mega_metric_matrix(71:72,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('69')

% N0
p.N0 = p.N0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.N0 = p.N0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.N0 = p.N0/(1-factor);%revert to old parameter

mega_metric_matrix(73:74,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('71')

% T0
p.T0 = p.T0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.T0 = p.T0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.T0 = p.T0/(1-factor);%revert to old parameter

mega_metric_matrix(75:76,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('73')
%
% p_L_I
p.p_L_I = p.p_L_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_I = p.p_L_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_I = p.p_L_I/(1-factor);%revert to old parameter

mega_metric_matrix2(1:2,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('1, second')

% eta_L_I
p.eta_L_I = p.eta_L_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_I = p.eta_L_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_I = p.eta_L_I/(1-factor);%revert to old parameter

mega_metric_matrix2(3:4,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('3, second')

% p_L_MPhi
p.p_L_MPhi = p.p_L_MPhi*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_MPhi = p.p_L_MPhi/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_MPhi = p.p_L_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix2(5:6,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('5, second')

% eta_L_MPhi
p.eta_L_MPhi = p.eta_L_MPhi*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_MPhi = p.eta_L_MPhi/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_MPhi = p.eta_L_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix2(7:8,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('7, second')

% p_L_M
p.p_L_M = p.p_L_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_M = p.p_L_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_L_M = p.p_L_M/(1-factor);%revert to old parameter

mega_metric_matrix2(9:10,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('9, second')

% eta_L_M
p.eta_L_M = p.eta_L_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_M = p.eta_L_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_L_M = p.eta_L_M/(1-factor);%revert to old parameter

mega_metric_matrix2(11:12,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('11, second')

% k_lin_L
p.k_lin_L = p.k_lin_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_L = p.k_lin_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_L = p.k_lin_L/(1-factor);%revert to old parameter

mega_metric_matrix2(13:14,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('13, second')

% k_B_L
p.k_B_L = p.k_B_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_L = p.k_B_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_L = p.k_B_L/(1-factor);%revert to old parameter

mega_metric_matrix2(15:16,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('15, second')

% k_U_L
p.k_U_L = p.k_U_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_L = p.k_U_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_L = p.k_U_L/(1-factor);%revert to old parameter

mega_metric_matrix2(17:18,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('17, second')

% k_int_L
p.k_int_L = p.k_int_L*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_L = p.k_int_L/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_L = p.k_int_L/(1-factor);%revert to old parameter

mega_metric_matrix2(19:20,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('19, second')

% p_G_MPhi_I
p.p_G_MPhi_I = p.p_G_MPhi_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_G_MPhi_I = p.p_G_MPhi_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_G_MPhi_I = p.p_G_MPhi_I/(1-factor);%revert to old parameter

mega_metric_matrix2(21:22,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('21, second')

% eta_G_MPhi
p.eta_G_MPhi = p.eta_G_MPhi*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_G_MPhi = p.eta_G_MPhi/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_G_MPhi = p.eta_G_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix2(23:24,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('23, second')

% p_G_M
p.p_G_M = p.p_G_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_G_M = p.p_G_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_G_M = p.p_G_M/(1-factor);%revert to old parameter

mega_metric_matrix2(25:26,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('25, second')

% eta_G_M
p.eta_G_M = p.eta_G_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_G_M = p.eta_G_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_G_M = p.eta_G_M/(1-factor);%revert to old parameter

mega_metric_matrix2(27:28,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('27, second')

% k_lin_G
p.k_lin_G = p.k_lin_G*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_G = p.k_lin_G/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_G = p.k_lin_G/(1-factor);%revert to old parameter

mega_metric_matrix2(29:30,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('29, second')

% k_B_G
p.k_B_G = p.k_B_G*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_G = p.k_B_G/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_G = p.k_B_G/(1-factor);%revert to old parameter

mega_metric_matrix2(31:32,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('31, second')

% k_U_G
p.k_U_G = p.k_U_G*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_G = p.k_U_G/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_G = p.k_U_G/(1-factor);%revert to old parameter

mega_metric_matrix2(33:34,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('33, second')

% k_int_G
p.k_int_G = p.k_int_G*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_G = p.k_int_G/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_G = p.k_int_G/(1-factor);%revert to old parameter

mega_metric_matrix2(35:36,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('35, second')

% p_C_M
p.p_C_M = p.p_C_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_C_M = p.p_C_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_C_M = p.p_C_M/(1-factor);%revert to old parameter

mega_metric_matrix2(37:38,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('37, second')

% eta_C_M
p.eta_C_M = p.eta_C_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_C_M = p.eta_C_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_C_M = p.eta_C_M/(1-factor);%revert to old parameter

mega_metric_matrix2(39:40,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('39, second')

% k_lin_C
p.k_lin_C = p.k_lin_C*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_C = p.k_lin_C/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_C = p.k_lin_C/(1-factor);%revert to old parameter

mega_metric_matrix2(41:42,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('41, second')

% k_B_C
p.k_B_C = p.k_B_C*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_C = p.k_B_C/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_C = p.k_B_C/(1-factor);%revert to old parameter

mega_metric_matrix2(43:44,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('43, second')

% k_U_C
p.k_U_C = p.k_U_C*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_C = p.k_U_C/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_C = p.k_U_C/(1-factor);%revert to old parameter

mega_metric_matrix2(45:46,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('45, second')

% k_int_C
p.k_int_C = p.k_int_C*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_C = p.k_int_C/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_C = p.k_int_C/(1-factor);%revert to old parameter

mega_metric_matrix2(47:48,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('47, second')

% p_F_I
p.p_F_I = p.p_F_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_I = p.p_F_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_I = p.p_F_I/(1-factor);%revert to old parameter

mega_metric_matrix2(49:50,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('49, second')

%
% eta_F_I
p.eta_F_I = p.eta_F_I*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_I = p.eta_F_I/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_I = p.eta_F_I/(1-factor);%revert to old parameter

mega_metric_matrix2(51:52,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('51, second')

% p_F_MpHI
p.p_F_MPhi = p.p_F_MPhi*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_MPhi = p.p_F_MPhi/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_MPhi = p.p_F_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix2(53:54,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('53, second')

% eta_F_MPhi
p.eta_F_MPhi = p.eta_F_MPhi*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_MPhi = p.eta_F_MPhi/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_MPhi = p.eta_F_MPhi/(1-factor);%revert to old parameter

mega_metric_matrix2(55:56,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('55, second')

% p_F_M
p.p_F_M = p.p_F_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_M = p.p_F_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.p_F_M = p.p_F_M/(1-factor);%revert to old parameter

mega_metric_matrix2(57:58,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('57, second')
%
% eta_F_M
p.eta_F_M = p.eta_F_M*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_M = p.eta_F_M/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.eta_F_M = p.eta_F_M/(1-factor);%revert to old parameter

mega_metric_matrix2(59:60,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('59, second')

% k_lin_F
p.k_lin_F = p.k_lin_F*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_F = p.k_lin_F/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_lin_F = p.k_lin_F/(1-factor);%revert to old parameter

mega_metric_matrix2(61:62,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('61, second')

% k_B_F
p.k_B_F = p.k_B_F*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_F = p.k_B_F/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_B_F = p.k_B_F/(1-factor);%revert to old parameter

mega_metric_matrix2(63:64,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('63, second')

% k_U_F
p.k_U_F = p.k_U_F*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_F = p.k_U_F/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_U_F = p.k_U_F/(1-factor);%revert to old parameter

mega_metric_matrix2(65:66,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('65, second')

% k_int_F
p.k_int_F = p.k_int_F*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_F = p.k_int_F/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.k_int_F = p.k_int_F/(1-factor);%revert to old parameter

mega_metric_matrix2(67:68,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('67, second')

% L_U0
p.L_U_0 = p.L_U_0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.L_U_0 = p.L_U_0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.L_U_0 = p.L_U_0/(1-factor);%revert to old parameter

mega_metric_matrix2(69:70,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('69, second')

% L_U0
p.G_U_0 = p.G_U_0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.G_U_0 = p.G_U_0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.G_U_0 = p.G_U_0/(1-factor);%revert to old parameter

mega_metric_matrix2(71:72,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('71, second')

% C_U0
p.C_U_0 = p.C_U_0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.C_U_0 = p.C_U_0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.C_U_0 = p.C_U_0/(1-factor);%revert to old parameter

mega_metric_matrix2(73:74,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('73, second')

% F_U0
p.F_U_0 = p.F_U_0*(1+factor); % increase by 10%
metric_vectorPLMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.F_U_0 = p.F_U_0/(1+factor)*(1-factor);
metric_vectorPLMPhi_10minus = sensitivity_calculator(p,tol,timeORIG,solORIG);
p.F_U_0 = p.F_U_0/(1-factor);%revert to old parameter

mega_metric_matrix2(75:76,:) = [metric_vectorPLMPhi_10plus;metric_vectorPLMPhi_10minus];
disp('75, second')


%% plot the most sensitive parameters

figure
hold on
ub_vec = max(mega_metric_matrix);
lb_vec = min(mega_metric_matrix);
 

col = [165,0,38;...
215,48,39;...
244,109,67;...
253,174,97;...
254,224,144;...
255,255,191;...
224,243,248;...
171,217,233;...
116,173,209;...
69,117,180;...
49,54,149]/255;%jet(11);

col_days = [64,0,75;...
118,42,131;...
153,112,171;...
194,165,207;...
231,212,232;...
247,247,247;...
217,240,211;...
166,219,160;...
90,174,97;...
27,120,55;...
0,68,27]/255;

fig1 = figure
fig2 = figure
fig3 = figure
fig4 = figure
%p 
figure(fig1)
hold on
j = 1;
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(1,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(1,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(1,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(1,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
figure(fig2)
hold on
j=1;
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(2,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(2,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(2,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(2,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%del_V_N
j = 2;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(9,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(9,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(9,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(9,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=2;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(10,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(10,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(10,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(10,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%beta
j = 3;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(11,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(11,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(11,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(11,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=3;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(12,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(12,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(12,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(12,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%eps_F_I
j = 4;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(3,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(3,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(3,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(3,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=4;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(4,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(4,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(4,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(4,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%del_I_T
j = 5;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(23,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(23,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(23,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(23,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=5;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(24,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(24,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(24,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(24,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_MPhi_I_L
j = 6;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(39,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(39,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(39,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(39,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=6;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(40,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(40,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
 for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(40,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(40,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_M_I
j = 7;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(45,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(45,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(45,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(45,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=7;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(46,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(46,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(46,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(46,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_T_I
j = 8;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(59,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(59,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(59,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(59,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=8;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix(60,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(60,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix(60,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(60,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_L_MPhi
j = 9;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(5,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(5,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(5,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(5,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=9;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(6,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(6,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(6,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(6,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_B_L
j = 10;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(15,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(15,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(15,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(15,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=10;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(16,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(16,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(16,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(16,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_int_L
j = 11;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(19,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(19,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(19,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(19,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=11;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(20,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(20,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(20,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(20,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_B_C
j = 12;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(43,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(43,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(43,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(43,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=12;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(44,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(44,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(44,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(44,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_F_I
j = 13;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(49,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(49,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(49,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(49,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=13;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(50,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(50,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(50,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(50,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_F_MPhi
j =14;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(53,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(53,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(53,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(53,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=14;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(54,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(54,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(54,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(54,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_lin_F
j = 15;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(61,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(61,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(61,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(61,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=15;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(62,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(62,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(62,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(62,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end

%k_B_F
j = 16;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(63,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(63,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(63,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(63,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=16;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(64,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(64,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(64,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(64,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end


%k_int_F
j = 17;
figure(fig1)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(67,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(67,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig3)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(67,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(67,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=17;
figure(fig2)
hold on
for i = 1:8
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),11);
   if isnan(mega_metric_matrix2(68,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(68,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
figure(fig4)
hold on
for i = 9:10
    grid_metric_1 = linspace(min([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),max([mega_metric_matrix(:,i);mega_metric_matrix2(:,i)]),9);
   if isnan(mega_metric_matrix2(68,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(68,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end


%p, del_V_N beta delta I T p_MPhi_L epsilon_L_MPhi p_T_I p_L_MPhi
%eta_L_MPhi p_L_M eta_L_M k_B_L k_U_L k_int_L p_F_I p_F_MPhi k_lin_F k_B_F
%k_int_F

figure(fig1)
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''}) %positive
    
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
xtickangle(45)
set(gca,'FontSize',18)

ylim([0 16])
xlim([0 8])

figure(fig2)
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{'p','\delta_{V,N}','\beta','\epsilon_{F,I}',...
    '\delta_{I,T}','p_{M\Phi,L}','p_{M,I}','p_{T,I}','p_{L,M\Phi}',...
  'k_{B_L}','k_{int_L}','k_{B_C}','p_{F,I}',...
   'p_{F,M\Phi}','k_{lin_F}','k_{B_F}','k_{int_F}'}) % negative
    
set(gca,'Xtick',[linspace(0.5,7.5,8)],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure'})%
xtickangle(45)
set(gca,'FontSize',18)

ylim([0 16])
xlim([0 8])


figure(fig3)
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'8','9'})
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{''})
    
set(gca,'Xtick',[9,10],'Xticklabel',{'Duration of tissue under 30%','Peak IFN time'})%%
xtickangle(45)
set(gca,'FontSize',18)
ylim([0 16])

figure(fig4)
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'8','9'})
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{'p','\delta_{V,N}','\beta','\epsilon_{F,I}',...
    '\delta_{I,T}','p_{M\Phi,L}','p_{M,I}','p_{T,I}','p_{L,M\Phi}',...
  'k_{B_L}','k_{int_L}','k_{B_C}','p_{F,I}',...
   'p_{F,M\Phi}','k_{lin_F}','k_{B_F}','k_{int_F}'})
    
set(gca,'Xtick',[9,10],'Xticklabel',{,'Duration of tissue under 30%','Peak IFN time'})%
xtickangle(45)
set(gca,'FontSize',18)
ylim([0 16])


figure
colormap(col);
c =  colorbar;
c.Ticks = [0 0.5 1];
c.TickLabels = {'max decrease','no change','max increase'};
set(gca,'FontSize',18)

figure
colormap(col_days);
c =  colorbar;
c.Ticks = [0 0.5 1];
c.TickLabels = {'max decrease','no change','max increase'};
set(gca,'FontSize',18)

STOP
full_sensitivity_analaysis_plot(mega_metric_matrix,mega_metric_matrix2)


 %%
 
%% plot the most sensitive parameters

figure
hold on
ub_vec = max(mega_metric_matrix);
lb_vec = min(mega_metric_matrix);
 

col = jet(11);
col_days = [255,247,243;...
253,224,221;...
252,197,192;...
250,159,181;...
247,104,161;...
221,52,151;...
174,1,126;...
122,1,119;...
73,0,106]/255; 

%p
j = 1;
for i = 1:8
   grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(1,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(1,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
   grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(1,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(1,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=2;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(2,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(2,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(2,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(2,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%del_V_N
j = 3;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(9,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(9,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(9,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(9,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=4;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(10,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(10,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(10,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(10,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%beta
j = 5;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(11,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(11,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(11,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(11,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=6;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(12,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(12,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(12,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(12,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%eps_F_I
j = 7;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(3,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(3,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(3,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(3,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=8;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(4,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(4,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(4,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(4,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%del_I_T
j = 9;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(23,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(23,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(23,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(23,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=10;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(24,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(24,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(24,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(24,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_MPhi_I_L
j = 11;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(39,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(39,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(39,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(39,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=12;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(40,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(40,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
 for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(40,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(40,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_M_I
j = 13;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(45,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(45,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(45,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(45,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=14;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(46,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(46,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(46,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(46,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_T_I
j = 15;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(59,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(59,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(59,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(59,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=16;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
   if isnan(mega_metric_matrix(60,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(60,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
   if isnan(mega_metric_matrix(60,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(60,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_L_MPhi
j = 17;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(5,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(5,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(5,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(5,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=18;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(6,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(6,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(6,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(6,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_B_L
j = 19;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(15,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(15,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(15,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(15,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=20;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(16,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(16,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(16,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(16,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_int_L
j = 21;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(19,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(19,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(19,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(19,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=22;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(20,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(20,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(20,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(20,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_B_C
j = 23;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(43,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(43,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(43,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(43,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=24;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(44,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(44,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(44,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(44,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_F_I
j = 25;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(49,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(49,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(49,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(49,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=26;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(50,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(50,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(50,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(50,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%p_F_MPhi
j = 27;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(53,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(53,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(53,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(53,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=28;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(54,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(54,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(54,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(54,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
%k_lin_F
j = 29;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(61,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(61,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(61,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(61,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=30;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(62,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(62,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(62,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(62,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end

%k_B_F
j = 31;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(63,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(63,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(63,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(63,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=32;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(64,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(64,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(64,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(64,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end


%k_int_F
j = 33;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(67,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(67,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(67,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(67,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end
j=34;
for i = 1:8
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
   if isnan(mega_metric_matrix2(68,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(68,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end
for i = 9:10
    grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
   if isnan(mega_metric_matrix2(68,i))==1
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix2(68,i),1);
       fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
   end
end


%p, del_V_N beta delta I T p_MPhi_L epsilon_L_MPhi p_T_I p_L_MPhi
%eta_L_MPhi p_L_M eta_L_M k_B_L k_U_L k_int_L p_F_I p_F_MPhi k_lin_F k_B_F
%k_int_F


set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,33.5,34),'Yticklabel',{'p^+','p^-','\delta_{V,N}^+','\delta_{V,N}^+','\beta^+','\beta^-','\epsilon_{F,I}^+','\epsilon_{F,I}^+',...
    '\delta_{I,T}^+','\delta_{I,T}^-','p_{M\Phi,L}^+','p_{M\Phi,L}^-','p_{M,I}^+','p_{M,I}^-','p_{T,I}^+','p_{T,I}^-','p_{L,M\Phi}^+','p_{L,M\Phi}^-',...
   'k_{B_L}^+','k_{B_L}^-', 'k_{int_L}^+','k_{int_L}^-','k_{B_C}^+','k_{B_C}^-','p_{F,I}^+','p_{F,I}^-',...
    'p_{F,M\Phi}^+','p_{F,M\Phi}^-','k_{lin_F}^+','k_{lin_F}^-','k_{B_F}^+','k_{B_F}^-','k_{int_F}^+','k_{int_F}^-'})
    
set(gca,'Xtick',[linspace(0.5,7.5,8),9,10],'Xticklabel',{'Max virus','Max dead cells','Minimum tissue','Max inflam. macs','Max T cells','Max IL-6','Max IFN','IFN exposure','Duration of tissue under 20%','Peak IFN time'})%
xtickangle(45)
set(gca,'FontSize',18)

ylim([0 34])
xlim([0 10.5])

    
STOP
 %%
col = jet(11);
col_days = [255,247,243;...
253,224,221;...
252,197,192;...
250,159,181;...
247,104,161;...
221,52,151;...
174,1,126;...
122,1,119;...
73,0,106]/255; 

% First figure subset
figure
hold on
ub_vec = max(mega_metric_matrix);
lb_vec = min(mega_metric_matrix);
for j = 1:24
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
       end
    end
    for i =9:10
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
       end
    end
end
for j = 41:41+7
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1 i-1 i i],[j-1-16 j-16 j-16 j-1-16],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1 i-1 i i],[j-1-16 j-16 j-16 j-1-16],col(val,:),'EdgeColor','none') 
       end
    end
    for i =9:10
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-16 j-16 j-16 j-1-16],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-16 j-16 j-16 j-1-16],col_days(val,:),'EdgeColor','none') 
       end
    end
end
set(gca,'Xtick',[linspace(0.5,8.5,9),10,11],'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,size(mega_metric_matrix,1)-0.5,size(mega_metric_matrix,1)),'Yticklabel',{'p^+','p^-','\epsilon_{F,I}^+','\epsilon_{F,I}^-','d_V^+','d_V^-',...
    '\delta_{V,M\Phi}^+','\delta_{V,M\Phi}^+','\delta_{V,N}^+','\delta _{V,N}^-','\beta^+','\beta^-','\lambda_S^+','\lambda_S^-','\delta_N^+','\delta_N^-',...
   'IC_{50,N}^+','IC_{50,N}^-','d_I^+','d_I^-','\delta_{I,M\Phi}^+','\delta_{I,M\Phi}^-','\delta_{I,T}^+','\delta_{I,T}^-','\epsilon_{G,M}^+','\epsilon_{G,M}^-',...
    '\psi_M^{max+}','\psi_M^{max-}','p_{M,I}^+','p_{M,I}^-','\epsilon_{I,M}^+','\epsilon_{I,M}^+'});%,...
   % 'V0^+','V0^-','S0^+','S0^-','M_{\Phi R,0}^+','M_{\Phi R,0}^-','M0^+','M0^-','N0^+','N0^-','T0^+','T0^-'})

set(gca,'Xtick',[linspace(0.5,8.5,8),9,10],'Xticklabel',{'max(V)','max(D)','min(S+R)','max(M_{\Phi I})','max(T)','max(L_U)','max(F_U)','F_U exposure','time ((S+R)/S_{max}<0.2)','day max(F_U)'})%
xtickangle(60)
set(gca,'FontSize',15)
xlim([0 10.5])
ylim([0 32])


col = jet(11);
figure
hold on
ub_vec = max(mega_metric_matrix);
lb_vec = min(mega_metric_matrix);
for j = 25:40%size(mega_metric_matrix,1)
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1 i-1 i i],[j-1-24 j-24 j-24 j-1-24],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1 i-1 i i],[j-1-24 j-24 j-24 j-1-24],col(val,:),'EdgeColor','none') 
       end
    end
   for i = 9:10
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-24 j-24 j-24 j-1-24],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-24 j-24 j-24 j-1-24],col_days(val,:),'EdgeColor','none') 
       end
    end
end
for j = 49:size(mega_metric_matrix,1)-12
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),11);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1 i-1 i i],[j-1-32 j-32 j-32 j-1-32],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1 i-1 i i],[j-1-32 j-32 j-32 j-1-32],col(val,:),'EdgeColor','none') 
       end
    end
    for i = 9:10
        grid_metric_1 = linspace(min(mega_metric_matrix(:,i)),max(mega_metric_matrix(:,i)),9);
       if isnan(mega_metric_matrix(j,i))==1
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-32 j-32 j-32 j-1-32],'k','EdgeColor','none') 
       else
           val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
           fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-32 j-32 j-32 j-1-32],col_days(val,:),'EdgeColor','none') 
       end
    end
end
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,size(mega_metric_matrix,1)-0.5,size(mega_metric_matrix,1)),'Yticklabel',{'\delta_{M\Phi,D}^+','\delta_{M\Phi,D}^-','\delta_{D,M\Phi}^+','\delta_{D,M\Phi}^-','a_{I,M\Phi}^+','a_{I,M\Phi}^-',...
    '\epsilon_{V,M\Phi}^+','\epsilon_{V,M\Phi}^-','\epsilon_{G,M\Phi}^+','\epsilon_{G,M\Phi}^-','p_{M\Phi,G}^+','p_{M\Phi,G}^-'...
    '\epsilon_{L,M\Phi}^+','\epsilon_{L,M\Phi}^-','p_{M\Phi,L}^+','p_{M\Phi,L}^-','\psi_N^{max+}','\psi_N^{max-}'...
    '\epsilon_{C,N}^+','\epsilon_{C,N}^-','p_{N,L}^+','p_{N,L}^-','\epsilon_{L,N}^+','\epsilon_{L,N}^-',...
    '\epsilon_{L,T}^+','\epsilon_{L,T}^-','p_{T,I}^+','p_{T,I}^-','\epsilon_{F,T}^+','\epsilon_{F,T}^-','p_{T,F}^+','p_{T,F}^-'});
set(gca,'Xtick',[linspace(0.5,8.5,8),9,10],'Xticklabel',{'max(V)','max(D)','min(S+R)','max(M_{\Phi I})','max(T)','max(L_U)','max(F_U)','F_U exposure','time ((S+R)/S_{max}<0.2)','day max(F_U)'})%
xtickangle(60)
set(gca,'FontSize',15)
ylim([0 40-24+size(mega_metric_matrix,1)-12-49+1])
xlim([0 10.5])

    figure
hold on
ub_vec = max(mega_metric_matrix2);
lb_vec = min(mega_metric_matrix2);
for j = 1:36%size(mega_metric_matrix2,1)
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
        if isnan(mega_metric_matrix2(j,i))==1
            fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
        else
            val = find(grid_metric_1>=mega_metric_matrix2(j,i),1);
            fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
        end
    end
        for i = 9:10
        grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
        if isnan(mega_metric_matrix2(j,i))==1
            fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],'k','EdgeColor','none') 
        else
            val = find(grid_metric_1>=mega_metric_matrix2(j,i),1);
            fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1 j j j-1],col_days(val,:),'EdgeColor','none') 
        end
    end
end
set(gca,'Xtick',linspace(0.5,10.5,11),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,size(mega_metric_matrix2,1)-0.5,size(mega_metric_matrix2,1)),'Yticklabel',{'p_{L,I}^+','p_{L,I}^-',...
    '\eta_{L,I}^+','\eta_{L,I}^-','p_{L,_M\Phi}^+','p_{L,M\Phi}^-','\eta_{L,M\Phi}^+','\eta_{L,M\Phi}^-','p_{L,_M}^+','p_{L,M}^-','\eta_{L,M}^+','\eta_{L,M}^-',...
    'k_{lin}_L^+','k_{lin}_L^-','k_B_L^+','k_B_L^+','k_U_L^+','k_U_L^-','k_{int}_L^+','k_{int}_L^-','p_{G,M\Phi}^+','p_{GL,M\Phi}^-','\eta_{G,M\Phi}^+','\eta_{G,M\Phi}^-',...
    'p_{G,M}^+','p_{G,M}^-','\eta_{G,M}^+','\eta_{G,M}^-','k_{lin}_G^+','k_{lin}_G^-','k_B_G^+','k_B_G^+','k_U_G^+','k_U_G^-','k_{int}_G^+','k_{int}_G^-'});%,...
   % 'L_{U,0}^+','L_{U,0}^-','G_{U,0}^+','G_{U,0}^-','C_{U,0}^+','C_{U,0}^-','F_{U,0}^+','F_{U,0}^-'})
set(gca,'Xtick',[linspace(0.5,8.5,8),9,10],'Xticklabel',{'max(V)','max(D)','min(S+R)','max(M_{\Phi I})','max(T)','max(L_U)','max(F_U)','F_U exposure','time ((S+R)/S_{max}<0.2)','day max(F_U)'})%
xtickangle(60)
set(gca,'FontSize',15)
ylim([0 36])
xlim([0 10.5])

figure
hold on
ub_vec = max(mega_metric_matrix2);
lb_vec = min(mega_metric_matrix2);
for j = 37:size(mega_metric_matrix2,1)
    for i = 1:8
        grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),11);
        if isnan(mega_metric_matrix2(j,i))==1
            fill([i-1 i-1 i i],[j-1-36 j-36 j-36 j-1-36],'k','EdgeColor','none') 
        else
            val = find(grid_metric_1>=mega_metric_matrix2(j,i),1);
            fill([i-1 i-1 i i],[j-1-36 j-36 j-36 j-1-36],col(val,:),'EdgeColor','none') 
        end
    end
     for i = 9:10
        grid_metric_1 = linspace(min(mega_metric_matrix2(:,i)),max(mega_metric_matrix2(:,i)),9);
        if isnan(mega_metric_matrix2(j,i))==1
            fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-36 j-36 j-36 j-1-36],'k','EdgeColor','none') 
        else
            val = find(grid_metric_1>=mega_metric_matrix2(j,i),1);
            fill([i-1+0.5 i-1+0.5 i+0.5 i+0.5],[j-1-36 j-36 j-36 j-1-36],col_days(val,:),'EdgeColor','none') 
        end
    end
end
set(gca,'Xtick',linspace(0.5,9.5,10),'Xticklabel',{'1','2','3','4','5','6','7','8','9','10','11'})
set(gca,'Ytick',linspace(0.5,size(mega_metric_matrix2,1)-0.5,size(mega_metric_matrix2,1)),'Yticklabel',{'p_{C,M}^+','p_{C,M}^-','\eta_{C,M}^+','\eta_{C,M}^-','k_{lin}_C^+','k_{lin}_C^-','k_B_C^+','k_B_C^+','k_U_C^+','k_U_C^-','k_{int}_C^+','k_{int}_C^-',...
    'p_{F,I}^+','p_{F,I}^-','\eta_{F,I}^+','\eta_{F,I}^-','p_{F,MPhi}^+','p_{F,MPhi}^-','\eta_{F,MPhi}^+','\eta_{F,MPhi}^-',...
    'p_{F,M}^+','p_{F,M}^-','\eta_{F,M}^+','\eta_{F,M}^-','k_{lin}_F^+','k_{lin}_F^-','k_B_F^+','k_B_F^+','k_U_F^+','k_U_F^-','k_{int}_F^+','k_{int}_F^-'});%,...
   % 'L_{U,0}^+','L_{U,0}^-','G_{U,0}^+','G_{U,0}^-','C_{U,0}^+','C_{U,0}^-','F_{U,0}^+','F_{U,0}^-'})

set(gca,'Xtick',[linspace(0.5,7.5,8),9,10],'Xticklabel',{'max(V)','max(D)','min(S+R)','max(M_{\Phi I})','max(T)','max(L_U)','max(F_U)','F_U exposure','time ((S+R)/S_{max}<0.2)','day max(F_U)'})%
xtickangle(60)
set(gca,'FontSize',15)
ylim([0 32])
xlim([0 10.5])

figure
colormap(col);
c =  colorbar;
c.Ticks = [0 0.5 1];
c.TickLabels = {'max neg','zero','max pos'};
set(gca,'FontSize',15)

figure
colormap(col_days);
c =  colorbar;
c.Ticks = [0 1];
c.TickLabels = {'min','max'};
set(gca,'FontSize',15)
