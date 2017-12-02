%% prepare for the computing environments
% author: Pengcheng Zhou, Columbia University, 2017
% email: zhoupc1988@gmail.com

%% get the directories.
ind = strfind(work_dir, 'code');
home_folder = work_dir(1:(ind-1));
code_folder = [home_folder, 'code', filesep]; 
cnmfe_folder = [home_folder, 'code', filesep, 'CNMF_E', filesep];
fig_folder = [home_folder, 'Figs', filesep];
video_folder = [home_folder, 'Videos', filesep];
data_folder = [home_folder, 'data', filesep];
results_folder = [home_folder, 'results', filesep];

%% configure CNMF-E
try
    temp = Sources2D();
    clear temp;
catch
    run([cnmfe_folder, 'cnmfe_setup.m']);
end
cd(cnmfe_folder);
[~, hash_cnmfe] = system('git rev-parse HEAD');
cd(work_dir);

%% set the seed of randomnuss generation for repeatable research
rng(1);