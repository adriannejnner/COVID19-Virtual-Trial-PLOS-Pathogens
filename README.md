# COVID19-Virtual-Trial-PLOS-Pathogens
This code accompanies the paper in PLOS Pathogens titled "COVID-19 virtual patient cohort suggests immune https://doi.org/10.1371/journal.ppat.1009753 mechanisms driving disease outcomes" 

To simulate the code that creates Figure 2 in the paper run "Commands_humandata_SIVD.m" in folder "CODE Figure 2 (SIVD refit)"

To simulate the code that creates Figure 3 onwards, you'll first need to install ddesd_f5, follow the instructions below: 
First download the dde solver: "ddesdinstallationv2.zip". To install the solver run "installMyDdesd()" in administrator mode in Matlab. Don't save this folder where you save the code

To simulate the code that creates Figure 3 in the paper run "commands_virus_IFN_model.m" in folder "CODE Figure 3 (IFN delay)"

To simulate the code that creates Figure 4 in the paper run "commands_model_full_and_disease_free_V2.m" in folder "CODE Figure 4 (Full model simulation)"

To simulate the code that creates Figure 5 in the paper run "Commands_sensitivity_analysis_V2.m" to create A and "Commands_sensitivity_analysis_QUALITATIVE_V2" for B-E in folder "CODE Figure 5 (Sensitivity analysis)"

To simulate the code that creates Figure 6 in the paper run "commands_beta_varying.m" in the folder "CODE Figure 6 (Beta change)"

To simulate the code that creates Figure 7 in the paper run "commands_fullmodel_macschanges.m" to simulate macrophage knockout, "commands_fullmodel_monochanges.m" to simulate monocyte knockout and "commands_fullmodel_neutschanges.m" to simulate neutrophil knockout in folder "CODE Figure 7 (knockout) and combine resulting figures. 

