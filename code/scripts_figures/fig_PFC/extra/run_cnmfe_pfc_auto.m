%% clear workspace
load /data/zhoupc/Data/paper_results/pfc_results.mat
neuron_manual = neuron.copy(); 
Coor_manual = Coor_cnmfe; 
neuron_manual_ds = neuron_ds.copy(); 
neuron = neuron_init.copy(); 
%% select data and map it to the RAM
cnmfe_choose_data;

%% create Source2D class object for storing results and parameters
Fs = 15;             % frame rate
ssub = 2;           % spatial downsampling factor
tsub = 2;           % temporal downsampling factor
gSig = 3;           % width of the gaussian kernel, which can approximates the average neuron shape
gSiz = 12;          % maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron_full.Fs = Fs;         % frame rate

% with dendrites or not 
with_dendrites = false;
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = gSiz;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end

%% options for running deconvolution 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold 

%% downsample data for fast and better initialization
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= 9000;             % user input: how many frames to read   (optional, default: until the end)

tic;
cnmfe_load_data;
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix

%%
neuron = neuron_init.copy(); 
neuron.options.min_pixel = 5; 
%% estimate the background 
% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1;  % spatial range / diameter of neurons

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
Nspatial = 10;       % this variable has different meanings:
%1) udpate_spatial_method=='hals' or 'hals_thresh',
%then Nspatial is the maximum iteration
%2) update_spatial_method== 'nnls', it is the maximum
%number of neurons overlapping at one pixel

neuron.options.maxIter = 2;   % iterations to update C

% parameters for running iteratiosn
nC = size(neuron.C, 1);    % number of neurons
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%%

%% pick neurons from the residual
neuron.options.min_corr = 0.65; % use a higher threshold for picking residual neurons.
neuron.options.min_pnr = 10;
seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
Yres = Ysignal - neuron.A*neuron.C; 
[center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, seed_method); % method can be either 'auto' or 'manual'
neuron_res = neuron.copy(); 

%% update background 
tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

%% update spatial and temporal 
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);


%% merge neurons 
merge_thr = [0.01, 0.6, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges

dmin = 1; 
center_method = 'max'; 
cnmfe_merge_neighbors; 
%% update spatial & temporal components
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center

    cnmfe_quick_merge;              % run neuron merges
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% order neurons 
snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2); 
% snr = max(neuron.C_raw, [], 2)./std(neuron.C_raw, 0, 2); 
[snr, srt] = sort(snr, 'descend'); 
neuron.orderROIs(srt); 
% neuron.viewNeurons([], neuron.C_raw); 

%% pick neurons from the residual
neuron.options.min_corr = 0.6; % use a higher threshold for picking residual neurons.
neuron.options.min_pnr = 10;
patch_par = [3,3];
seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
% [center_new, Cn_res, pnr_res] = 
neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, seed_method); % method can be either 'auto' or 'manual'

%% update background 
tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

%% update spatial and temporal 
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% merge neurons 
merge_thr = [0.01, 0.6, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges
% 
dmin = 0.5; 
center_method = 'max'; 
cnmfe_merge_neighbors; 
%% update spatial and temporal 
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
%% merge neurons 
merge_thr = [0.01, 0.6, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges



merge_thr = [0.8, 0.01, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges
% 
dmin = 1; 
center_method = 'max'; 
cnmfe_merge_neighbors; 

%% update spatial and temporal 
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
%% merge neurons 
merge_thr = [0.01, 0.6, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges



merge_thr = [0.8, 0.01, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges
% 
dmin = 1; 
center_method = 'max'; 
cnmfe_merge_neighbors; 

%% order neurons 
snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2); 
[snr, srt] = sort(snr, 'descend'); 
neuron.orderROIs(srt); 
neuron.viewNeurons([], neuron.C_raw); 

%% update background 
tic;
spatial_ds_factor = 1; 
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

%% update spatial and temporal 
tic;
for m=1:1
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
        %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    %         neuron.A = (Ysignal*neuron.C')/(neuron.C*neuron.C');
            neuron.trimSpatial(0.02, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    
    if isempty(merged_ROI)
        break;
    end
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% merge neurons 
merge_thr = [0.01, 0.6, 0.];
% merge neurons
cnmfe_quick_merge;              % run neuron merges
% 
dmin = 0.5; 
center_method = 'max'; 
cnmfe_merge_neighbors; 


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

