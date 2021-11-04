% simulating different beta values

% model parameters initialised

tspan = [0 30];

% simulate model for original 0%
[time0,sol0] = COVID_IMMUNE_MODEL(p,tspan);
disp('0%')

% simulate model for 10% increase
p.beta =0.3*1.1;
[time10,sol10] = COVID_IMMUNE_MODEL(p,tspan);
disp('10%')

% simulate model for 20% increase
p.beta = 0.3*1.2;
[time20,sol20] = COVID_IMMUNE_MODEL(p,tspan);
disp('20%')

% simulate model for 30% increase
p.beta = 0.3*1.3;
[time30,sol30] = COVID_IMMUNE_MODEL(p,tspan);
disp('30%')

% simulate model for 40% increase
p.beta = 0.3*1.4;
[time40,sol40] = COVID_IMMUNE_MODEL(p,tspan);
disp('40%')

% simulate model for 50% increase
p.beta = 0.3*1.5;
[time50,sol50] = COVID_IMMUNE_MODEL(p,tspan);
disp('50%')

beta_plotted(time0,sol0,time10,sol10,time20,sol20,time30,sol30,time40,sol40,time50,sol50,p)

