%% simulated annealing to create new virtual cohort


%% Quickly checking how many old patients do actually satisfy new viral load bounds?
% from old patients, do we have to exclude any? if so how many

% load old patients


% simulate their dynamics


% are they within the bounds, yes or no?
rng default % for reproducivibility


%% Creating new patients


%p = load_parameters_simulated_annealing2();

%number of patients
N = 100;

%calculating beta std
std_day3_Munster = 1.1958;
mean_day3_Munster = 5.9966;
relative_std_day3_Munster = std_day3_Munster/mean_day3_Munster;
sigma_beta = relative_std_day3_Munster; 

%calculating eps_F_I std
sigma_eps_F_I = 208.9/(3125000); %standard deviation sigma informed by confidence intervals for data fit

% calculating p_F_I,eta_F_MPhi and eta_F_I std
CI_day0_TrouilletAssant = [2*1e2 3*1e3]/1000; %pg/ml
mean_day0_TrouilletAssant = 3.5*1e2/1000;
number_patients_TrouilletAssant = 26; 
t_score_26df_95CI = 2.060;
std_TrouilletAssant = (CI_day0_TrouilletAssant(2)-mean_day0_TrouilletAssant)/t_score_26df_95CI*sqrt(number_patients_TrouilletAssant);
relative_std_day0_TrouilletAssant = std_TrouilletAssant/mean_day0_TrouilletAssant;
sigma_p_F_I = relative_std_day0_TrouilletAssant/10;
sigma_eta_F_MPhi = relative_std_day0_TrouilletAssant/1e6*5;
sigma_eta_F_I = relative_std_day0_TrouilletAssant/1e3;

%calculating p_L_MPhi and p_MPhi_L std
mean_mild_Herold = 21;
std_mild_Herold = 19;
mean_severe_Herold = 195;
std_severe_Herold = 165;

mean_moderate_Lucas = 2.449560839602535e+02;
std_moderate_Lucas = 5.385157964298745e+02;
relative_std_moderate_Lucas = std_moderate_Lucas/mean_moderate_Lucas;

mean_severe_Lucas = 1.004274975130175e+03;
std_severe_Lucas = 1.646578521236555e+03;
relative_std_severe_Lucas = std_severe_Lucas/mean_severe_Lucas;

relative_std_mild_Herold = std_mild_Herold/mean_mild_Herold;
relative_std_severe_Herold = std_severe_Herold/mean_severe_Herold;

sigma_p_L_MPhi_mild = relative_std_moderate_Lucas;%relative_std_mild_Herold;
sigma_p_MPhi_I_L_mild = relative_std_moderate_Lucas;%relative_std_mild_Herold*p.p_MPhi_I_L/10;
sigma_p_L_MPhi_severe = relative_std_severe_Lucas;%relative_std_severe_Herold;
sigma_p_MPhi_I_L_severe = relative_std_severe_Lucas;%relative_std_severe_Herold*p.p_MPhi_I_L/10;

%calculating eps_L_T and tau_T std
mean_Liu_Tcells_day10 = 0.35;
std_Liu_Tcells_day10 = 0.1;
relative_std_Liu_Tcells_day10 = std_Liu_Tcells_day10/mean_Liu_Tcells_day10;
sigma_eps_L_T = relative_std_Liu_Tcells_day10*p.eps_L_T;
sigma_tau_T = relative_std_Liu_Tcells_day10;

% calculating p_M_I std
mean_Liu_monocytes_day4 = 1;
std_Liu_monocytes_day4 = 0.75;
relative_std_Liu_monocytes_day4 = std_Liu_monocytes_day4/mean_Liu_monocytes_day4;
sigma_p_M_I = relative_std_Liu_monocytes_day4*p.p_M_I/2;

%calculating p_F_M std
sigma_p_F_M = 126.2*1e-4; %standard deviation sigma informed by confidence intervals for data fit

mu_vec =    [p.beta,    p.p_MPhi_I_L,    p.p_L_MPhi,    p.p_F_I,    p.eta_F_I,    p.eps_L_T,    p.p_M_I,    p.eta_F_MPhi,    p.tau_T     p.eps_F_I,    p.p_F_M];
sigma_vec = [sigma_beta,sigma_p_MPhi_I_L_mild,sigma_p_L_MPhi_mild,sigma_p_F_I,sigma_eta_F_I,sigma_eps_L_T,sigma_p_M_I,sigma_eta_F_MPhi,sigma_tau_T,sigma_eps_F_I,sigma_p_F_M];

%     figure
% for i = 1:length(mu_vec)
%     subplot(3,4,i)
%     hold on 
%     hist(normrnd(mu_vec(i),sigma_vec(i),1,200))
% end
% 
% 
% lb = max(repmat(0,1,length(mu_vec)),mu_vec-5*sigma_vec);
% ub = mu_vec+5*sigma_vec;
% 
% for i = 1:N
%     
%      x0 = normrnd(mu_vec,sigma_vec);%initial guess vector
%         
%      if isempty(find(x0<0))==0 
%         while isempty(find(x0<0))==0
%             x0 = normrnd(mu_vec,sigma_vec/2);%initial guess vector
%         end
%      end
%      
%      options = optimoptions('simulannealbnd','ObjectiveLimit',0.1);%','PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf},
%  
%     [x,fval,exitFlag,output] = simulannealbnd(@(x)COVID_minimising_fun_mod(x,p),x0,lb,ub,options);
%     patient(i,:) = [x];
%     disp('patient number')
%     i
% end


mu_vec =    [p.beta,    p.p_MPhi_I_L,    p.p_L_MPhi,    p.p_F_I,    p.eta_F_I,    p.eps_L_T,    p.p_M_I,    p.eta_F_MPhi,    p.tau_T     p.eps_F_I,    p.p_F_M];
sigma_vec = [sigma_beta,sigma_p_MPhi_I_L_severe,sigma_p_L_MPhi_severe,sigma_p_F_I,sigma_eta_F_I,sigma_eps_L_T,sigma_p_M_I,sigma_eta_F_MPhi,sigma_tau_T,sigma_eps_F_I,sigma_p_F_M];
% 
%     figure
% for i = 1:length(mu_vec)
%     subplot(3,4,i)
%     hold on 
%     hist(normrnd(mu_vec(i),sigma_vec(i),1,200))
% end


lb = max(repmat(0,1,length(mu_vec)),mu_vec-5*sigma_vec);%mu_vec-5*sigma_vec);
ub = mu_vec+5*sigma_vec;%5*sigma_vec;

STOP

for i =1:4*N
    
     x0 = normrnd(mu_vec,sigma_vec);%initial guess vector
        
     if isempty(find(x0<0))==0 
        while isempty(find(x0<0))==0
            x0 = normrnd(mu_vec,sigma_vec/2);%initial guess vector
        end
     end
     
     options = optimoptions('simulannealbnd','ObjectiveLimit',0.1);%','PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf},
 
    [x,fval,exitFlag,output] = simulannealbnd(@(x)COVID_minimising_fun_severe(x,p),x0,lb,ub,options);
    patient(i,:) = [x];
    
 
    disp('patient number')
    i
end

 save('patients_sim_annealing_severe15.mat','patient')

 %%
 
 load('patients_sim_annealing_severeNEW2.mat')
figure
hold on 
hist(patient(:,2))
plot([mu_vec(2) mu_vec(2)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('p_{MPhi_I,L}')
set(gca,'FontSize',18)

figure
hold on 
hist(patient(:,3))
plot([mu_vec(3) mu_vec(3)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('p_{L,MPhi}')
set(gca,'FontSize',18)

figure
hold on 
hist(patient(:,4))
plot([mu_vec(4) mu_vec(4)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{F,I}')

figure
hold on 
hist(patient(:,7))
plot([mu_vec(7) mu_vec(7)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{M,I}')

figure
hold on 
hist(patient(:,8))
plot([mu_vec(8) mu_vec(8)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('\eta_{F,MPhi}')

figure
hold on 
hist(patient(:,10))
plot([mu_vec(10) mu_vec(10)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('\epsilon_{F,I}')

figure
hold on 
hist(patient(:,11))
plot([mu_vec(11) mu_vec(11)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{F,M}')

figure
hold on 
hist(patient(:,2))
plot([mu_vec(2) mu_vec(2)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('p_{MPhi_I,L}')
set(gca,'FontSize',18)

figure
hold on 
hist(patient(:,3))
plot([mu_vec(3) mu_vec(3)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('p_{L,MPhi}')
set(gca,'FontSize',18)

figure
hold on 
hist(patient(:,4))
plot([mu_vec(4) mu_vec(4)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{F,I}')

figure
hold on 
hist(patient(:,7))
plot([mu_vec(7) mu_vec(7)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{M,I}')

figure
hold on 
hist(patient(:,8))
plot([mu_vec(8) mu_vec(8)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('\eta_{F,MPhi}')

figure
hold on 
hist(patient(:,10))
plot([mu_vec(10) mu_vec(10)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('\epsilon_{F,I}')

figure
hold on 
hist(patient(:,11))
plot([mu_vec(11) mu_vec(11)],[0 40],'--','Color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'FontSize',18)
xlabel('p_{F,M}')