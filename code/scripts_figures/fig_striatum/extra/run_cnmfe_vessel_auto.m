load /data/zhoupc/Data/paper_results/vessel_results.mat
neuron_manual = neuron.copy(); 
Coor_manual = Coor_cnmfe; 
neuron = neuron_init.copy(); 


%% order neurons 
snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2); 
[~, srt] = sort(snr, 'descend'); 
neuron.orderROIs(srt); 

%% estimate the background 
% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = true;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 30;       % this variable has different meanings:
%1) udpate_spatial_method=='hals' or 'hals_thresh',
%then Nspatial is the maximum iteration
%2) update_spatial_method== 'nnls', it is the maximum
%number of neurons overlapping at one pixel

neuron.options.maxIter = 2;   % iterations to update C
%% update spatial & temporal components
tic;
for m=1:1
    %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% merge neurons 

merge_thr = [0.01, 0.7, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges


dmin = 2; 
center_method = 'max'; 
cnmfe_merge_neighbors; 

%% pick neurons from the residual
neuron.options.min_corr = 0.8; % use a higher threshold for picking residual neurons.
neuron.options.min_pnr = 10;
patch_par = [3,3];
seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
[center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, seed_method); % method can be either 'auto' or 'manual'

%% update spatial & temporal components
tic;
for m=1:1
    %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% estimate the background 
% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

%% merge neurons 

merge_thr = [0.01, 0.7, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges

%% update spatial & temporal components
tic;
for m=1:1
    %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% order neurons 
snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2); 
[snr, srt] = sort(snr, 'descend'); 
neuron.orderROIs(srt); 
neuron.viewNeurons([], neuron.C_raw); 

%% pick neurons from the residual
seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
Yres = Ysignal - neuron.A*neuron.C; 
[center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, seed_method); % method can be either 'auto' or 'manual'
neuron_res = neuron.copy(); 


%% apply results to the full resolution
if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    spatial_ds_factor = tsub * spatial_ds_factor; 
    cnmfe_full;
    neuron_full = neuron.copy();
end

%% order neurons 
snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2); 
[snr, srt] = sort(snr, 'descend'); 
neuron.orderROIs(srt); 
neuron.viewNeurons([], neuron.C_raw); 

%% display neurons
dir_neurons = sprintf('%s%s%s_neurons%s', dir_nm, filesep, file_nm, filesep);
if exist('dir_neurons', 'dir')
    temp = cd();
    cd(dir_neurons);
    delete *;
    cd(temp);
else
    mkdir(dir_neurons);
end
neuron.viewNeurons([], neuron.C_raw, dir_neurons);
close(gcf); 

%% display contours of the neurons
figure;
% Cn_nobg = correlation_image(neuron.reshape(Ysignal(:, 1:5:end), 2), 8);
neuron.Coor = plot_contours(neuron.A, Cn_nobg, 0.6, 0, [], [], 1);
colormap gray;
axis equal; axis off;
title('contours of estimated neurons');

% plot contours with IDs
% [Cn_filter, pnr] = neuron.correlation_pnr(Ysignal);
figure;
plot_contours(neuron.A, Cn_filter, 0.6, 0, [], neuron.Coor, 1);
colormap gray;
title('contours of estimated neurons');

%% check spatial and temporal components by playing movies
% save_avi = false;
% avi_name = 'play_movie.avi';
% neuron.Cn = Cn;
% neuron.runMovie(Ysignal, [0, 50], save_avi, avi_name);
kt = 1;     % play one frame in every kt frames
t_begin = 1; %neuron.Fs*60*3; 
t_end = T; %neuron.Fs*60*5; 
save_avi = true;
y_quantile = 0.9999;    % for specifying the color value limits 
ac_quantile = 0.9995;

cnmfe_save_video;

%% visually check the results by choosing few neurons
neuron.runMovie(Ysignal, [0, max(Ysignal(:))]*0.8); 
%% save results
globalVars = who('global');
file_nm(file_nm==' ') = '_'; 
eval(sprintf('save %s%s%s_results.mat %s', dir_nm, filesep, file_nm, strjoin(globalVars)));

