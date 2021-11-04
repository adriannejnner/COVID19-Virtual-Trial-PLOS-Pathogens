function [max_IL6 max_neutrophils minimum_tissue, max_inflam_macs, IFN_exposure, IFN_peak, max_Tcells, max_monocytes]=reevaluate_patients(time_grid,sol_virus,sol_tissue,sol_infected,sol_dead,sol_macs_res,sol_macs_inflam,sol_neutrophils, sol_monocytes,sol_Tcells, sol_IL6, sol_IFN, sol_GCSF, time)



max_IL6 = max(sol_IL6');
max_neutrophils = max(sol_neutrophils');
max_Tcells = max(sol_Tcells');
max_inflam_macs = max(sol_macs_inflam');
max_monocytes = max(sol_monocytes');

minimum_tissue = min(sol_tissue');

tgrid_vec = [];
for i = 1:999
   tgrid_vec(i) = time(i+1)-time(i); 
end

for ii = 1:size(sol_virus,1)
    peak_loc = find(sol_IFN(ii,:)== max(sol_IFN(ii,:)));
    IFN_peak(ii) = time(peak_loc); %time of IFN peak
    IFN_exposure(ii) = sum(sol_IFN(ii,1:end-1).*tgrid_vec);
end


end