function p = load_parameters()

% LIST OF PARAMETERS

p.rho = 0.5;   % bystander damage modulation

% VIRUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.V0         = 12;         % initial virus Copies/mL
p.eps_F_I    = 133;%625;             % IFN inhibition of viral production
p.del_V_MPhi = 768;       % Rate macrophages clear viruses (1/day)
p.del_V_N    = p.del_V_MPhi;    % Rate neutrophils clear viruses (1/day)
p.d_V        = 18.94;               % decay rate of viruses (turned off when immune presence is considered)
p.phat       = 741.2;            %lytic viral production rate

% CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.R0 = 0;

%Epithelial cell %------------------------------------
p.S0        = 0.16;         % initial alveolar type I +alveolar type II + epithelial cells x 10^9
p.lam_S     = 0.7397;           % epithelial cell replication ate
p.Smax	    = p.S0;             % carrying capacity
p.h_N       = 3.0198;           % Hill coefficient for neutrophil induced damage of cells (susceptible and infected)
p.IC_50_N   = 4.7054*1e-5;           % IC50 for neutrophil induced damage of cells (susceptible and infected)
p.del_N     = 1.6786; 	        % Rate neutrophils inflict tissue damage (1/cell/day)- XX needs to be converted to ng/mL/day (susceptible and infected)
p.beta      = 0.289;        %Infection rate
p.tau_I     = 0.1667;             % converted - eclipse phase

%Infected cell %------------------------------------
p.I0         = 0;               % Initial infected cells x 10^9
p.del_I_MPhi = 121.195;    % Rate macrophages phagocytose infected cells (1/day)
p.del_I_T    = 50000;            % Rate T cells clear infected cell (1/cell 1/day)
p.d_I	     = 0.144;            % lysis rate

%Dead/damaged cell %------------------------------------
p.D0         = 0;               % initial dead cells x 10^9
p.del_D_MPhi = 8.0256;      %rate macrophages phagocytose dead cells (1/cell 1/day)
p.d_D        = 8;               % converted - dead cell degradation rate (1/day)

%Resident macrophage %------------------------------------
p.MPhi_R_0     = 23*1e-3/843;             % Initial inactive macrophages x 10^8 cells XXXX THIS USED TO BE 36!!
p.d_MPhi_R     = 0.0124;	     % death rate of resident macrophages
p.lam_MPhi     = 5.9432*1e3;     % production or resident macrophages from infiltrating (systemic) macs
p.eps_V_MPhi   = 905.2205;       % half-effect virus concentration for resident macrophage replenishing
p.MPhi_max     = p.MPhi_R_0;     % carrying capacity for resident macs
p.a_I_MPhi     = 1.1*1e3;         % activation of resident macrophages to inflammatory macros
p.del_MPhi_D   = 6.0607;     % rate macs die from phaocytosis

%Infiltrating macrophage %------------------------------------
p.MPhi_I_0     = 0;             % Initial infiltrating macrophages x 10^8 cells 
p.p_MPhi_I_G   = 0.42;        % differentiation of monocyte to macrophage by GM-CSF (rate)
p.eps_G_MPhi   = 2664.5;        % differentation of monocyte to macrophage by GM-CSf (half-effect)
p.h_M_MPhi     = 2.0347;     	% differentation of monocyte to macrophage by GM-CSf (hill coefficient)
p.p_MPhi_I_L   = p.p_MPhi_I_G;             % differentiation of monocyte to macrophage by IL-6 (rate)
p.eps_L_MPhi   = 1102.9;        % differentiation of monocyte to macrophage by IL-6 (half-effect)
p.d_MPhi_I     = 0.3;        % death rate of infiltrating macrophages

%Monocytes %--------------------------------------------------------------
p.M0           = 0.0004;                % Initial monocytes x 10^8
p.eps_G_M      = 57.197;            % GM-CSF recruitement of monocytes (half-effect)
p.d_M          = 0.7562;            % death rate of monocytes (1/day)
p.h_M          = 1.6711;
p.p_M_I        = 2.2*1e-1;       % recruitment of monocytes by infected cells                            UNKNOWN                                                                                                      
p.eps_I_M      = p.Smax/3;          % half effect recruitment of monocytes                                  UNKNOWN     
p.MR           = 0.1619*70/5000; 
p.psi_M_max    = 11.5451;

% Neutrophils %-----------------------------------------------------------
p.N0         = 0.00526; %!         % Initial neutrophils x 10^8
p.NR         = 0.03156; %!    % Number of neutrophils in the resevoir x 10^8
p.eps_C_N    = 1.8924*1e-4;     % G-CSF recruitment of neutrophils (half-effect)
p.d_N        = 1.2773;%0.8;     % Death rate of neutrophils (1/day) 
p.psi_N_max  = 4.1335;          % maximal neutrophil production
p.p_N_L      = 0.3;           % production of neutrophils by IL-6 (rate)                              UNKNOWN!!!
p.eps_L_N    = p.eps_G_M;       % production of neutrophils by IL-6 (half-effect)

% T cells %-----------------------------------------------------------
p.T0         = 1.104*1e-4;           % Initial number ot T cells x()??
p.p_T_L      = 4;           %production of T cells by IL-6 (production rate)
p.eps_L_T    = 474.49;      % production of T cells by IL-6 (half-effect)
p.p_T_F      = 4;           %production of T cells by IFN (production rate)
p.eps_F_T    = 399.3;       % production of T cells by IFN (half-effect)
p.d_T        = 0.4;         % decay rate of T cells
p.p_T_I      = 1;           % recruitment of T cells by infected cells  (rate)
p.eps_T_I    = 1e-6;         % recruitment of T cells by infected cells (half-effect) 
p.tau_T      = 4.5;         % delay in T cell arrival

% CYTOKINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specific cytokine parameters
p.stoch= 1;                 % stochiometric constant
p.stoch_C = 1.4608;         % stochiometric constant for G-CSF
p.avo=6.02214E23;           % avogadro's number

% IL-6 --------------------------------------------------------
p.L_U_0 = 1.1;             % Initial GM-CSF (unbound) in ng/mL
p.L_B_0 = 0;                % Initial GM-CSF (bound) in ng/mL
p.p_L_I = 1188.7;
p.eta_L_I = 7.232*1e-9;
p.p_L_MPhi = 1872;          %IL-6 production by infiltrating (systemic) macrophages (rate)
p.eta_L_MPhi = 10^-5;        %IL-6 production by infiltrating (systemic) macrophages (half-effect)
p.p_L_M = 7.256*1e4;               %IL-6 production by monocytes (rate)
p.eta_L_M = 8.99*1e4/1e9;               %IL-6 production by monocytes (half-effect)
p.k_lin_L = 16.636;%14.6872;        % renal clearance
p.k_int_L = 61.8;            % internalisation rate
p.k_B_L = 0.0018;           % (converted) binding rate
p.k_U_L = 22.29;           % (converted) unbinding rate
p.R_L_N = 720;             % Number of IL-6 receptors on neutrophils
p.R_L_T = 300;             % Number of IL-6 receptors on T cells
p.R_L_M = 509;             % Number of IL-6 receptors on monocytes
p.MM_L = 21000;              % Molecular weight of GM-CSF

% GM-CSF --------------------------------------------------------
p.G_U_0 = 2.43;             % Initial GM-CSF (unbound) in ng/mL
p.G_B_0 = 0;                % Initial GM-CSF (bound) in ng/mL
p.p_G_MPhi_I = 20000;         % production of GM-CSF by macrophages                                 UNKNOWN
p.eta_G_MPhi = 0.1;   %unknown   % half-effect production GM-CSF by monocytes                      UNKNOWN BUT ASSUMED EQUAL TO P_G_M AT HOMEOSTASIS
p.p_G_M = 10680;              % production of GM-CSF by monocytes
p.eta_G_M = 0.14931;   %unknown   % half-effect production GM-CSF by monocytes                         UNKNOWN BUT ESTIMATED FROM HOMEOSTASIS
p.k_lin_G = 11.7427;        % renal clearance
p.k_int_G = 73.44;            % internalisation rate
p.k_B_G = 0.0021;           % (converted) binding rate
p.k_U_G = 522.72;           % (converted) unbinding rate
p.R_G_M = 1058;%1820;             % Number of GM-CSF receptors on monocytes   % XXXXXXXXX CHANGES THIS FROM 1820 as couldn't find original reference
p.MM_G = 14E3;              % Molecular weight of GM-CSF

% G-CSF --------------------------------------------------------
p.C_U_0 = 0.025;       % Initial G-CSF (unbound) in ng/mL
p.C_B_0 = 2.1899E-5;   % Initial G-CSF (bound) in ng/mL
p.p_C_M = p.p_G_M/1000; %unknown             % production rate of G-CSF by monocyte XXXXXXXXXXXXXX                           UNKNOWN
p.eta_C_M = p.M0*2;     %unknown          % half-effect production rate of G-CSF by monocytes  XXXXXXXXXXXXXXXXXX         UNKNOWN - ESTIMATE FROM HOMEOSTASIS?
p.k_lin_C = 0.16139;        % renal clearance
p.k_int_C = 462.42;         % internalisation rate
p.k_B_C = 2.243;           % binding rate (ng/mL)^-1day^-1
p.k_U_C = 184.87;           % (converted) unbinding rate
p.R_C_N = 600;              % Number of G-CSF receptors on neutrophils
p.MM_C = 19.6E3;            % Molecular weight of G-CSF

% IFN -----------------------------------------------------------
p.F_U_0 = 1.5*1e-2;
p.F_B_0 = 0.307482;
p.p_F_I = 2.823*1e4;
p.eta_F_I = 0.001116;
p.p_F_MPhi = 1;                                                                                    %UNKNOWN
p.eta_F_MPhi = p.MPhi_R_0*2;        %                                                              %UNKNOWN
p.p_F_M = 3.56*1e4;
p.eta_F_M = 5.405*1e4/1e9;
p.k_lin_F = 16.635;
p.k_int_F = 16.968;
p.k_B_F = 0.0107;
p.k_U_F = 6.072;
p.R_F_T = 1000;
p.R_F_I = 1300;
p.MM_F = 19000;

%----- BINDING CONVERSION FACTOR --------------------------------------------

p.A_L = p.stoch*p.MM_L/p.avo*(p.R_L_M+p.R_L_N+p.R_L_T)*(1/5000)*10^9*1e12;% IL-6
p.A_G = p.stoch*p.MM_G/p.avo*(p.R_G_M)*(1/5000)*10^9*1e12;% GM-CSF
p.A_C = 2*p.MM_C/p.avo*p.R_C_N*(1/5000)*10^9*1e9;% G-CSF
p.A_F = p.stoch*p.MM_F/p.avo*(p.R_F_T+p.R_F_I)*(1/5000)*10^9*1e12;% IFN


%% EXTRA

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
%p.p_G_M = 7.7*1e5;
%p.eps_G_MPhi = 1.846806147120571e-05;%p.eps_G_MPhi/(p.G_U_0/p.G_B_0)

%modulated IFN parameters to get desired IFN dynamics in the IFN submodel
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;

p.p_F_M = 3.5600*1e4*1e-4; %changed

%fixing IL-6 range
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;%1e-2*1.3;%1102.9/1e7;% HERE!!!!!-
p.eps_G_MPhi = 2664.5/1e5;%1e-3*1.3;%2664.5/1e7;%%1.846806147120571e-05;

p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;%0.5*1e-3/0.5;

p.p_T_I = 0.01*0.8; %changed from 1
p.del_I_T = 238;%50;%50;238 got from T cell killing paper

p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
%
%p.tau_T = 6;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; % CHANGED FROM 5
p.eps_I_M = 0.11;
p.del_V_N = 768*3;%*1.5;
p.del_V_MPhi = 768*100;%10000*4*10;%*1.5;

%p.d_MPhi_I = 1;
%p.k_lin_C = 16;%

%p.d_D = 21;
%p.d_M = p.d_N;
%p.p_MPhi_I_G = 0.42*2;
%p.p_MPhi_I_L = 0.42*2;

p.p_G_M = 1.234*1e3*1e-1;
%p.p_G_MPhi_I = 1232.4;

p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%p.eta_G_MPhi = 1e-3;

% XXXXXXXXXXXXXXX
% p.d_V = 18.94*1.2;
% p.d_V_spec = 8;
% p.d_I = 0.0144*6;
% p.phat = 741.2*2.5;
% p.beta = 0.289*0.6;
% p.del_V_N = 768*3*2;%*1.5;
% p.del_V_MPhi = 768*100*2;%10000*4*10;%*1.5;
% XXXXXXXXXXXXXXX
%%

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
%p.p_G_M = 7.7*1e5;
%p.eps_G_MPhi = 1.846806147120571e-05;%p.eps_G_MPhi/(p.G_U_0/p.G_B_0)

%modulated IFN parameters to get desired IFN dynamics in the IFN submodel
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;

p.p_F_M = 3.5600*1e4*1e-4; %changed

%fixing IL-6 range
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;%1e-2*1.3;%1102.9/1e7;% HERE!!!!!-
p.eps_G_MPhi = 2664.5/1e5;%1e-3*1.3;%2664.5/1e7;%%1.846806147120571e-05;

p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;%0.5*1e-3/0.5;

p.p_T_I = 0.01*0.8; %changed from 1
p.del_I_T = 238;%50;%50;238 got from T cell killing paper

p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
%
%p.tau_T = 6;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; % CHANGED FROM 5
p.eps_I_M = 0.11;
p.del_V_N = 768*3;%*1.5;
p.del_V_MPhi = 768*100;%10000*4*10;%*1.5;

%p.d_MPhi_I = 1;
%p.k_lin_C = 16;%

%p.d_D = 21;
%p.d_M = p.d_N;
%p.p_MPhi_I_G = 0.42*2;
%p.p_MPhi_I_L = 0.42*2;

p.p_G_M = 1.234*1e3*1e-1;
%p.p_G_MPhi_I = 1232.4;

p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%p.eta_G_MPhi = 1e-3;

% XXXXXXXXXXXXXXX
% p.d_V = 18.94*1.2;
% p.d_V_spec = 8;
% p.d_I = 0.0144*6;
% p.phat = 741.2*2.5;
% p.beta = 0.289*0.6;
% p.del_V_N = 768*3*2;%*1.5;
% p.del_V_MPhi = 768*100*2;%10000*4*10;%*1.5;
% XXXXXXXXXXXXXXX

p = Homeostasis_calculations(p);

end