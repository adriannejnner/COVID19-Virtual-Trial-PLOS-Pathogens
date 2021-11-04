function y = COVID_minimising_fun_severe(x,p)

   % p = load_parameters_simulated_annealing;
    
   
    p.p_MPhi_I_L = x(2);
    p.p_L_MPhi = x(3);
    p.p_F_I = x(4);
    p.p_M_I = x(7);
    p.eta_F_MPhi =  x(8);
    p.eps_F_I =  x(10);
    p.p_F_M = x(11);
    
    p = Homeostasis_calculations(p);
    estimated_params = [p.F_B_0,p.G_B_0,p.C_B_0,p.L_B_0,p.MPhi_I_0,p.eta_L_MPhi,p.p_C_M,p.eta_C_M,p.p_N_L];

    if isempty(find(estimated_params<0))==0
        disp('Negative parameter')
        y = 1e3;
    elseif isempty(find(estimated_params>1e9))==0
        disp('Extremely large parameter')
        y = 1e3;
    else
        
        tspan = [0 15];

        if isempty(find(x<0))==0 
           disp('One of the parameters is negative-BAD') 
        end

        % evaluate model
        [time,solSS,sol] = COVID_IMMUNE_MODELINSILICO(p,tspan);

        if isreal(sol.y)==0
            y = 1e4;
        else     
        %evaluate objective function    
            load('output_ranges.mat')
            load('IL6_Lucas_data.mat')
            load('viral_load_boundsNEW.mat')%','lower_St','upper_St')
            
            ub_IFN = ub_IFN*0.1;
            lb_IFN(1) = 0.1;
            lb_IFN = lb_IFN*0.1;

            ub_IL6 = [400 400 400];%[200 200 200];%max(10.^IL6_Lucas_Severe_individual);%40;
            lb_IL6 = [0 0 0];
            time_IL6 = [4.5 7 10];
            lower_St(1) = 2;
            
            vec_virus = sum(max((deval(sol,1:15,1)-(lower_St(1:15)+upper_St(1:15))/2).^2-(upper_St(1:15)-(lower_St(1:15)+upper_St(1:15))/2).^2,0)); %previous only used 1-7
            vec_IFN = sum(max((deval(sol,time_IFN,17)-(lb_IFN+ub_IFN)/2).^2-(ub_IFN-(lb_IFN+ub_IFN)/2).^2,0));
            %vec_IL6 = sum(max((deval(sol,time_IL6,11)-(lb_IL6+ub_IL6)/2).^2-(ub_IL6-(lb_IL6+ub_IL6)/2).^2,0));
            vec_IL6 = sum(max((deval(sol,time_IL6,11)-(lb_IL6+ub_IL6)/2).^2-(ub_IL6-(lb_IL6+ub_IL6)/2).^2,0));
            vec_GCSF = sum(max((deval(sol,time_GCSF,13)-(lb_GCSF+ub_GCSF)/2).^2-(ub_GCSF-(lb_GCSF+ub_GCSF)/2).^2,0));
            vec_neutrophils = 0;%sum(max((deval(sol,time_neutrophils,9)-(lb_neutrophils+ub_neutrophils)/2).^2-(ub_neutrophils-(lb_neutrophils+ub_neutrophils)/2).^2,0))*10;

             y = [vec_virus+vec_IFN+vec_IL6+vec_GCSF+vec_neutrophils]
        end
    end
end