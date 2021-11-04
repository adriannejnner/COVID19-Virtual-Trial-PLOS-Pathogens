function metric_vectorEpsVMPhi_10plus = sensitivity_calculator(p,tol,timeORIG,solORIG)

p = Homeostasis_calculations(p);% check model is still at homeostasis

tspan = [0 50];
[timeH, solH] = COVID_IMMUNE_MODEL_disease_free_full_model(p,tspan);
p.C_U_0 = p.C_U_0*1.2;
[timeP, solP] = COVID_IMMUNE_MODEL_disease_free_full_model(p,tspan);
p.C_U_0 = p.C_U_0/1.2;

if sum(abs(solP.y(:,1)-solP.y(:,end)))<tol %still at homeostasis
    tspan = [0 30];
    
    %XXXX Changed from COVID_IMMUNE_MODEL to COVID_IMMUNE_MODEL_FM
    [timeEpsVMPhi_10plus,solEpsVMPhi_10plus] = COVID_IMMUNE_MODEL_FM(p,tspan);% solve model
    metric_vectorEpsVMPhi_10plus = metrics(timeEpsVMPhi_10plus,solEpsVMPhi_10plus,timeORIG,solORIG);% calculate all the metrics
    if isreal(metric_vectorEpsVMPhi_10plus)==0
        metric_vectorEpsVMPhi_10plus = NaN(1,10);
    end
else
    disp('not at homeostasis')
    metric_vectorEpsVMPhi_10plus = NaN(1,10);
end

hold on 
plot(timeEpsVMPhi_10plus,solEpsVMPhi_10plus(17,:))
end