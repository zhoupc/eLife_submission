%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
cellsort_folder = [code_folder, 'CellSort'];
addpath(cellsort_folder);

output_folder = [fig_folder, filesep, 'Fig_SIM_subfigs'];
results_folder = sprintf('%sfig_corr%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_corr.mat'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 8, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 0.6, ...   % GB, space for loading data within one patch
    'patch_dims', [100, 100]);  %GB, patch size

use_parallel = true;
update_sn = true;
save_initialization = false;

%% load the background and noise levels extracted from the a typical microendoscopic data
load data_BG;
[d1,d2] = size(sn);
T0 = size(F,2);
Fs0 = 15;   % frame rate used in the background
D = reshape(D, d1*d2, []);
sn = reshape(sn,d1*d2, 1);

%% parameters for simulating neurons
K = 200;    % number of neurons
gSig = 3;   % gaussian width of cell body
gSiz = 4*gSig+1;
d_min = 2*gSig;  % minimum distance between neurons
mu_bar = 0.1;  % mean spike counts per frame
minIEI = 10;    % minimum interval between events, unit: frame
minEventSize = 2;  % minimum event size
T = 2000;   % number of frames
Fs = 10;     % frame rate
neuron_amp = 2.0;   % amplify neural signals to modify the signal to noise ratio
if Fs~=Fs0
    F = resample(F', Fs, Fs0)';
    T0 = size(F,2);
end
if T<T0
    F = F(:, 1:T);
else
    F = resample(F', T, T0)';
end
seed = 2;
tau_d_bar = 1*Fs;  % mean decay time constant, unit: frame
tau_r_bar = 0.2*Fs;  % mean rising time constant, unit: frame

%% simulate cellular signals
sim_AC_corr;
ind_center = sub2ind([d1,d2], round(coor.y0), round(coor.x0));
A = bsxfun(@times, A, reshape(sn(ind_center), 1, []));   %adjust the absolute amplitude neuron signals based on its spatial location

%% select two neurons only
cell_id2 = [139, 152];
A_2 = A(:, cell_id2);

temp = cos_similarity(C');
[idx12, ~, ~] = find(temp==min(abs(temp(:))));
[~, idx3] = min(sum(temp(idx12,:)));
C_2 = C([idx12; idx3], :);
S_2 = S([idx12; idx3], :);
img = reshape(sum(A_2, 2), d1, d2);
[r, c] = find(img>0);
rmin = min(r)-30;
rmax = max(r)+30;
cmin = min(c)-30;
cmax = max(c)+30;
ind = false(d1,d2);
d1_2 = rmax-rmin+1;
d2_2 = cmax-cmin+1;
ind(rmin:rmax, cmin:cmax) = true;
A2 = A_2(ind, :);
D2 = D(ind,:);
sn = sn(ind);

%% generate background and noise
E = bsxfun(@times, randn(d1_2*d2_2, T), sn);
Bc = D2*mean(F,2);       % time-invariant baseline
c3 = C_2(3,:);
c3 = (c3-mean(c3))*std(F(1,:))/std(c3);
F(1,:) = F(1,:)-c3;
Bf = D2*bsxfun(@minus, F, mean(F,2));  % fluctuating background

%% run CNMF-E given different correlation
neuron = Sources2D('d1', d1_2, 'd2', d2_2);
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
debug_on = false;

min_corr = 0.8;     % minimum local correlation for a seeding pixel
min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 30;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
bg_neuron_ratio = 1.5;
ring_radius = 18;

neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd, 'gSig', gSig, 'gSiz', gSiz,...
    'deconv_flag', true, 'ring_radius', 18, 'background_model', 'ring');
neuron.options.deconv_options = struct('type', 'ar2', ...
    'method', 'foopsi', ...
    'optimize_pars', true, ...
    'smin', -5);
neuron.Fs = Fs0;
neuron0 = neuron.copy();


%% find the values that generates the desired correlation
x = linspace(0, 1, 100);
y = zeros(size(x));
for m=1:length(x)
    c = x(m);
    temp = corr(C_2(1:2,:)'*(1-c) + C_2([3,3],:)'*c);
    y(m) = temp(2);
end
gamma_val = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
c_val = interp1(y, x, gamma_val);

%% try different correlation coefficients;
A_all = cell(length(gamma_val), 1);
C_all = A_all;
C_raw = C_all;
spatial_sim_all = A_all;
temporal_sim_all = C_all;
ind_match_all = A_all;
nPCs = 20;
nICs = 20;

for mm = length(gamma_val):-1:1
    %% simulate data
    fprintf('corr. coeff = %.2f\n', gamma_val(mm));
    c = c_val(mm);
    neuron = neuron0.copy();
    N = 5;
    C2 = C_2(1:2, :)*(1-c)+C_2([3,3], :) * c;
    A_est = zeros(d1_2*d2_2, 2, N);
    C_est = zeros(2, T, N);
    C_raw_est = zeros(2, T, N);
    spatial_sim = nan(2, N);
    temporal_sim = nan(2, N);
    ind_match = nan(2, N);
    Y = round(neuron_amp*A2*C2 + bsxfun(@plus, Bf, Bc) + E);  % observed fluorescence
    Y = neuron.reshape(Y, 2);
    Ysiz = size(Y);
    nam = sprintf('%sdata_%.2f.mat', results_folder, gamma_val(mm));
    save(nam, 'Y', 'Ysiz', '-v7.3');
    neuron.select_data(nam);
    
    %% run PCA/ICA analysis
    rng(1);
    [A_ica, C_ica] = run_pca_ica(neuron.reshape(Y, 2), nPCs, nICs);
    delete('temp.tif'); 
    rmdir('temp', 's'); 
    neuron.A = A_ica;
    neuron.C_raw = C_ica;
    neuron.C = C_ica;
    neuron.compactSpatial();
    match_results = pair_neurons(A2, C2, neuron.A, neuron.C);
    
    for m=1:2
        ids = match_results.ind_max(m);
        if isnan(ids)
            continue;
        end
        A_est(:,m, N) = neuron.A(:,ids);
        C_est(m,:,N) = neuron.C(ids, :);
        C_raw_est(m,:,N) = neuron.C_raw(ids, :);
        ind_match(m,N) = m;
        
        spatial_sim(m, N) = cos_similarity(A2(:, m), neuron.A(:, ids));
        temporal_sim(m, N) = cos_similarity(neuron.C(ids,:)', C2(m, :)');
    end
    
    %% run CNMF-E
    rng(1);
    neuron.getReady(pars_envs);
    
    [center, Cn, PNR] = neuron.initComponents_parallel(K, [], save_initialization, use_parallel);
    match_results = pair_neurons(neuron.A, neuron.C, A2, C2);
    kk = size(neuron.A,2);
    for m=1:kk
        ids = match_results.ind_spatial(m);
        if isnan(ids)
            continue;
        end
        A_est(:,m, 1) = neuron.A(:,m);
        C_est(m,:,1) = neuron.C(m, :);
        C_raw_est(m,:,1) = neuron.C_raw(m, :);
        ind_match(m,1) = ids;
        ind(ids) = false;
        spatial_sim(m, 1) = cos_similarity(A2(:, ids), neuron.A(:, m));
        temporal_sim(m, 1) = cos_similarity(neuron.C(m,:)', C2(ids, :)');
    end
    %% estimate the background components
    neuron.update_background_parallel(use_parallel);
    neuron_init = neuron.copy();
    
    %% pick neurons from the residual
    if kk<2
        neuron.initComponents_residual_parallel(2-kk, save_initialization, use_parallel, 0.5, 6);
    end
    match_results = pair_neurons(neuron.A, neuron.C, A2, C2);
    for m=1:2
        ids = match_results.ind_spatial(m);
        if isnan(ids)
            continue;
        end
        A_est(:,m, 2) = neuron.A(:,m);
        C_est(m,:,2) = neuron.C(m, :);
        C_raw_est(m,:,2) = neuron.C_raw(m, :);
        ind_match(m,2) = ids;
        
        spatial_sim(m, 2) = cos_similarity(A2(:, m), neuron.A(:, ids));
        temporal_sim(m, 2) = cos_similarity(neuron.C(ids,:)', C2(m, :)');
    end
    
    %% correct the inaccurate extraction of neural traces
    use_c_hat = false;
    neuron.update_temporal_parallel(use_parallel, use_c_hat);
    neuron.update_spatial_parallel(use_parallel);
    match_results = pair_neurons(neuron.A, neuron.C, A2, C2);
    for m=1:2
        ids = match_results.ind_spatial(m);
        if isnan(ids)
            continue;
        end
        A_est(:,m, 3) = neuron.A(:,m);
        C_est(m,:,3) = neuron.C(m, :);
        C_raw_est(m,:,3) = neuron.C_raw(m, :);
        ind_match(m,3) = ids;
        
        spatial_sim(m, 3) = cos_similarity(A2(:, m), neuron.A(:, ids));
        temporal_sim(m, 3) = cos_similarity(neuron.C(ids,:)', C2(m, :)');
    end
    %% udpate spatial&temporal components,
    for m=1:2
        % update temporal
        neuron.update_temporal_parallel(use_parallel);
        
        % update _spatial
        neuron.update_spatial_parallel(use_parallel);
    end
    match_results = pair_neurons(neuron.A, neuron.C, A2, C2);
    for m=1:2
        ids = match_results.ind_spatial(m);
        if isnan(ids)
            continue;
        end
        A_est(:,m, 4) = neuron.A(:,m);
        C_est(m,:,4) = neuron.C(m, :);
        C_raw_est(m,:,4) = neuron.C_raw(m, :);
        ind_match(m,4) = ids;
        
        spatial_sim(m, 4) = cos_similarity(A2(:, m), neuron.A(:, ids));
        temporal_sim(m, 4) = cos_similarity(neuron.C(ids,:)', C2(m, :)');
    end
    
    %     neuron.save_workspace();
    
    %% collect results
    A_all{mm} = A_est;
    C_all{mm} = C_est;
    C_raw_all{mm} = C_raw_est;
    spatial_sim_all{mm} = spatial_sim;
    temporal_sim_all{mm} = temporal_sim;
    ind_match_all{mm} = ind_match;
end
% end

%%  initialization results
a1_sim_init = nan(size(gamma_val));
a2_sim_init = a1_sim_init;
c1_sim_init = a1_sim_init;
c2_sim_init = a1_sim_init;

for mm=1:length(gamma_val)
    ind_match = ind_match_all{m};
    spatial_sim = spatial_sim_all{mm}(:, 2);
    temporal_sim = temporal_sim_all{mm}(:,2);
    for m=1:2
        idx = ind_match(m);
        if idx==1
            a1_sim_init(mm) = spatial_sim(m);
            c1_sim_init(mm) = temporal_sim(m);
        elseif idx==2
            a2_sim_init(mm) = spatial_sim(m);
            c2_sim_init(mm) = temporal_sim(m);
        else
            continue;
        end
    end
end


%% final results
a1_sim_final = nan(size(gamma_val));
a2_sim_final = a1_sim_init;
c1_sim_final = a1_sim_init;
c2_sim_final = a1_sim_init;

for mm=1:length(gamma_val)
    ind_match = ind_match_all{m};
    spatial_sim = spatial_sim_all{mm}(:, end-1);
    temporal_sim = temporal_sim_all{mm}(:, end-1);
    for m=1:2
        idx = ind_match(m);
        if idx==1
            a1_sim_final(mm) = spatial_sim(m);
            c1_sim_final(mm) = temporal_sim(m);
        elseif idx==2
            a2_sim_final(mm) = spatial_sim(m);
            c2_sim_final(mm) = temporal_sim(m);
        else
            continue;
        end
    end
end
%% PCA/ICA results
a1_sim_ica = nan(size(gamma_val));
a2_sim_ica = a1_sim_ica;
c1_sim_ica = a1_sim_ica;
c2_sim_ica = a1_sim_ica;

for mm=1:length(gamma_val)
    ind_match = ind_match_all{m};
    spatial_sim = spatial_sim_all{mm}(:, end);
    temporal_sim = temporal_sim_all{mm}(:, end);
    for m=1:2
        idx = ind_match(m);
        if idx==1
            a1_sim_ica(mm) = spatial_sim(m);
            c1_sim_ica(mm) = temporal_sim(m);
        elseif idx==2
            a2_sim_ica(mm) = spatial_sim(m);
            c2_sim_ica(mm) = temporal_sim(m);
        else
            continue;
        end
    end
end
%% results
x = 0:(length(gamma_val)-1);
y = atanh([a1_sim_init; a1_sim_final; a1_sim_ica]);
dx = 0.25;
col = [0.5, 0.5, 0.5; 1, 0.5, 0.3; 0.52, 0.16, 0.8] ;
alpha_v = 0.3;

w_left = 0.12;
w = 0.85;
h_top = 0.92;
h = 0.17;

figure('papersize', [8, 4]);
init_fig;
set(gcf, 'defaultaxesfontsize', 18);
axes('position', [0.12, -0.01, w, 1.033]);
hold on;
b = bar(x, 0*y', 'grouped');
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);
xlim([-dx*2, x(end)+2*dx]);
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);
ylim([1, 2]);
axis off;
temp = legend('CNMF-E initialization', 'CNMF-E final', 'PCA/ICA', 'orientation', 'horizental'); %, 'location', 'northwest');
set(temp, 'fontsize', 13);


axes('position', [w_left, h_top-h, w, h*0.9]);
hold on;
b = bar(x, y', 'grouped');
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);

xlim([-dx*2, x(end)+2*dx]);
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);
set(gca, 'ylim', [0, 3.5])
set(gca, 'ytick', atanh([0.2, 0.9, 0.99]))
set(gca, 'yticklabel',([0.2, 0.9, 0.99]))
ylabel('$a_1$', 'Interpreter', 'latex');
box off;

% legend('init.', 'iter. 1', 'with correction', 'orientation', 'horizental', 'fontsize', 13);

% draw barplot of c1
y = atanh([c1_sim_init; c1_sim_final; c1_sim_ica]);
axes('position', [w_left, h_top-2*h, w, h*0.9]);
hold on;
b = bar(x, y', 'grouped');
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);

xlim([-dx*2, x(end)+2*dx]);
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);
set(gca, 'ylim', [0, 3.5])
set(gca, 'ytick', atanh([0.2, 0.9, 0.99]))
set(gca, 'yticklabel', ([0.2, 0.9, 0.99]))
ylabel('$c_1$', 'Interpreter', 'latex');
box off;


% draw barplot of a2
y = atanh([a2_sim_init; a2_sim_final; a2_sim_ica]);
axes('position', [w_left, h_top-3*h, w, h*0.9]);
hold on; b = bar(x, y', 'grouped');
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);

xlim([-dx*2, x(end)+2*dx]);
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);
set(gca, 'ylim', [0, 3.5])
set(gca, 'ytick', atanh([0.2, 0.9, 0.99]))
set(gca, 'yticklabel', [0.2, 0.9, 0.99])
ylabel('$a_2$', 'Interpreter', 'latex');
box off;

% draw barplot of c2
y = atanh([c2_sim_init; c2_sim_final; c2_sim_ica]);
axes('position', [w_left, h_top-4*h, w, h*0.9]);
hold on; b = bar(x, y', 'grouped');
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);

xlim([-dx*2, x(end)+2*dx]);
set(gca, 'xtick', x);
set(gca, 'xticklabel', gamma_val);
set(gca, 'ylim', [0, 3.5])
set(gca, 'ytick', atanh([0.2, 0.9, 0.99]))
set(gca, 'yticklabel', ([0.2, 0.9, 0.99]))
ylabel('$c_2$', 'Interpreter', 'latex');
box off;
xlabel('corr$(c_1, c_2)$', 'Interpreter', 'latex');


if export_fig
    saveas(gcf, sprintf('%s/similarity_corr.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/similarity_corr.pdf', output_folder));
end
%%


%% show the evolution of CNMF-E
gamma = 0.9;
A_est = A_all{gamma_val==gamma};
C_est = C_all{gamma_val==gamma};
ind_match = ind_match_all{gamma_val==gamma};
C2 = C_2(1:2, :)*(1-c_val(gamma_val==gamma))+C_2([3,3], :) * c_val(end);
A2_sort = A2(:, ind_match(:, end-1));
C2_sort = C2(ind_match(:, end-1), :);

figure('papersize', [5, 3]);
init_fig;
w_left = 0.01;
w1 = 0.15;
w2 = 0.65;
h = 0.22;
h_top = 0.93;
tmp_xlim = d2_2/2+[-1,1]*gSiz;
tmp_ylim = d1_2/2+[-1,1]*gSiz;

% ground truth
neuron.A = A2_sort;
neuron.C = C2_sort;

img = zeros(d1_2, d2_2, 3);
img(:, :, 1) = neuron.reshape(neuron.A(:,1), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left, h_top-h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

img = zeros(d1_2, d2_2, 3);
img(:, :, 2) = neuron.reshape(neuron.A(:,2), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left+0.92*w1, h_top-h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

axes('position', [w_left+w1*2.0, h_top-h, w2, h*1.19]);
plot(neuron.C(1, :)/max(neuron.C(1,:))*1.5, 'g'); hold on;
plot(neuron.C(2,:)/max(neuron.C(2,:))*1.5+1, 'r');
ybg = mean(Bf,1);
plot(ybg/max(ybg)*1.5+2, 'b')
axis off;

ylim([0, 3.5]);
xlim([0, 1200]);
% initialization
neuron.A = A_est(:, :, 2);
neuron.C = C_est(:, :, 2);
img = zeros(d1_2, d2_2, 3);
img(:, :, 1) = neuron.reshape(neuron.A(:,1), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left, h_top-2*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;


img = zeros(d1_2, d2_2, 3);
img(:, :, 2) = neuron.reshape(neuron.A(:,2), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);
axes('position', [w_left+0.92*w1, h_top-2*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

axes('position', [w_left+w1*2.0, h_top-2*h, w2, h*0.85]);
plot(neuron.C(1, :)/max(neuron.C(1,:))*1.5, 'g'); hold on;
plot(neuron.C(2,:)/max(neuron.C(2,:))*1.5+1, 'r');
axis off;

ylim([0, 2.5]);
xlim([0, 1200]);

% update spatial & temporoal
neuron.A = A_est(:, :, end-1);
neuron.C = C_est(:, :, end-1);

img = zeros(d1_2, d2_2, 3);
img(:, :, 1) = neuron.reshape(neuron.A(:,1), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left, h_top-3*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

img = zeros(d1_2, d2_2, 3);
img(:, :, 2) = neuron.reshape(neuron.A(:,2), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left+0.92*w1, h_top-3*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

axes('position', [w_left+w1*2.0, h_top-3*h, w2, h*0.85]);
plot(neuron.C(1, :)/max(neuron.C(1,:))*1.5, 'g'); hold on;
plot(neuron.C(2,:)/max(neuron.C(2,:))*1.5+1, 'r');
axis off;

ylim([0, 2.5]);
xlim([0, 1200]);

% PCA/ICA
neuron.A = A_est(:, ind_match(:, end-1), end);
neuron.C = C_est(ind_match(:, end-1), :, end);

img = zeros(d1_2, d2_2, 3);
img(:, :, 1) = neuron.reshape(neuron.A(:,1), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left, h_top-4*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

img = zeros(d1_2, d2_2, 3);
img(:, :, 2) = neuron.reshape(neuron.A(:,2), 2);
img(:, :, 3) = sum(img, 3);
img = uint16(img/max(img(:))*65536);

axes('position', [w_left+0.92*w1, h_top-4*h, w1, h*0.85]);
imagesc(65536-img);
axis equal;
xlim(tmp_xlim);
ylim(tmp_ylim); box on; set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

axes('position', [w_left+w1*2.0, h_top-4.17*h, w2, h*1.02]);
plot(neuron.C(1, :)/max(neuron.C(1,:))*1.5, 'g'); hold on;
plot(neuron.C(2,:)/max(neuron.C(2,:))*1.5+1, 'r');
axis off;
ylim([-0.5, 2.5]);
xlim([0, 1200]);


if export_fig
    saveas(gcf, sprintf('%s/example_corr.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_corr.pdf', output_folder));
end




%%
save(results_file, 'A_all', 'C_all', 'C_raw_all', 'spatial_sim_all', 'temporal_sim_all', 'ind_match_all', 'A2', 'C2', 'gamma_val', '-append');
















