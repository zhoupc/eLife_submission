work_dir = fileparts(mfilename('fullpath')); 
cd(work_dir); 

%% Figure 1: microendoscopic data contain large background signals with rapid fluctuations due to multiple sources. 
cd(work_dir); 
run ./scripts_figures/fig_intro/main.m 

%% Figure 2: CNMF-E can accurately separate and recover the background fluctuations in simulated data 
cd(work_dir); 
run ./scripts_figures/fig_bg/main.m

%% Figure 3: CNMF-E accurately initializes individual neurons' spatial and temporal components in simulated data 
cd(work_dir); 
run ./scripts_figures/fig_initialization/main.m 

%% Figure 4: SIMULATION 
cd(work_dir); 
run ./scripts_figures/fig_sim/main.m 

%% Figure 5: correlation 
cd(work_dir); 
run ./scripts_figures/fig_correlation/main.m 

%% Figure 6: striatum 