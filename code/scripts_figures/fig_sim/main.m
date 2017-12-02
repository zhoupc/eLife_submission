%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_SIM_subfigs'];
results_folder = sprintf('%sfig_sim%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_sim.mat'];
cellsort_folder = [code_folder, 'CellSort'];
addpath(cellsort_folder);

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% computation environment
pars_envs = struct('memory_size_to_use', 8, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 1, ...   % GB, space for loading data within one patch
    'patch_dims', [64, 64]);  %GB, patch size

use_parallel = true;
update_sn = true;
save_initialization = false;
use_prev = false;   % use previous initialization results 

%% load the extracted background and noise from the data
load ../fig_initialization/data_BG;
[d1,d2] = size(sn);
T0 = size(F,2);
Fs0 = 15;   % frame rate used in the background
neuron = Sources2D('d1', d1, 'd2', d2);
D = neuron.reshape(D, 1);
sn = neuron.reshape(sn, 1);

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
    F = imresize(F, [size(F,1), round(size(F,2)*Fs/Fs0)]);
    T0 = size(F,2);
end
if T<T0
    F = F(:, 1:T);
else
    temp = repmat(F, [1, ceil(T/T0)]);
    F = temp(:, 1:T);
end
seed = 2;
tau_d_bar = 1*Fs;  % mean decay time constant, unit: frame
tau_r_bar = 0.2*Fs;  % mean rising time constant, unit: frame
neuron.updateParams('gSig', gSig, 'gSiz', gSiz);
%% simulate data
% cellular signals
sim_AC;
ind_center = sub2ind([d1,d2], round(coor.y0), round(coor.x0));
A = bsxfun(@times, A, reshape(sn(ind_center), 1, []));   %adjust the absolute amplitude neuron signals based on its spatial location
C = neuron_amp*C;

% simulate white noise
E = bsxfun(@times, randn(d1*d2, T), sn);

% background
Bc = D*mean(F,2);       % time-invariant baseline
Bf = D*bsxfun(@minus, F, mean(F,2));  % fluctuating background

%% SNR factor = 1
Y = A*C + bsxfun(@plus, Bf, Bc) + E;  % observed fluorescence
Y = round(Y); 
RSS_ground_truth = sum(E(:).^2); 
var_Y = sum(var(Y, 0, 2))*size(Y, 2); 

%% CNMF-E options
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
min_corr = 0.9;     % minimum local correlation for a seeding pixel
min_pnr = 15;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 1;            % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
ring_radius = 18;
bg_ssub = 1;   % downsample the background estimation 
dmin = 5;   % merge neurons that are close
merge_thr = 0.85;   % merge neurons that are highly correlated.
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd, 'gSig', gSig, 'gSiz', gSiz,...
    'deconv_flag', true, 'ring_radius', ring_radius, ...
    'background_model', 'ring', 'deconv_flag', true, ...
    'dmin', dmin, 'merge_thr', merge_thr, 'bg_ssub', bg_ssub, 'dist', 4, ...
    'dist', 4);
neuron.options.deconv_options = struct('type', 'ar2', ...
    'method', 'foopsi', ...
    'optimize_pars', true, ...
    'optimize_b', true, ...
    'smin', -3);
neuron.options.spatial_constraints.circular = true; 
neuron.Fs = Fs0;
neuron_basic = neuron.copy(); 

%% show how RSS changes over different operations. Figure 4B
try
    neuron_init = results_data.neuron0_init;
    neuron = results_data.neuron0;
    D = results_data.D;
    F = results_data.F;
    A = results_data.A;
    C = results_data.C;
    sn = results_data.sn;
    RSS_cnmfe = results_data.RSS_cnmfe;
    RSS_ids = results_data.RSS_ids;
    labels = results_data.labels;
    label_meanings = results_data.label_meanings;
catch
%     panel_rss_true; 
    panel_rss;
end

%% normalize RSS 
% RSS_cnmfe = RSS./RSS_ground_truth; 
% RSS_true = RSS_true./RSS_ground_truth; 
RSS_cnmfe_norm = RSS_cnmfe./RSS_cnmfe(end); 
% RSS_true = RSS_init_true/var_Y; 

%%
figure('papersize', [8,6]);
init_fig;
set(gcf, 'defaultAxesfontsize', 25);
axes('position', [0.22, 0.192, 0.77, 0.78]);
x = RSS_ids;
plot(x, RSS_cnmfe_norm, '-k', 'linewidth', 3); hold on;
h1 = plot(x(labels==1), RSS_cnmfe_norm(labels==1), '^', 'markersize', 18, 'markerfacecolor', 'b'); hold on;
h2 = plot(x(labels==2), RSS_cnmfe_norm(labels==2), 's', 'markersize', 18, 'markerfacecolor', 'm');
h3 = plot(x(labels==3), RSS_cnmfe_norm(labels==3), 'd', 'markersize', 18, 'markerfacecolor', 'c');
h4 = plot(x(labels==4), RSS_cnmfe_norm(labels==4), 'v', 'markersize', 18, 'markerfacecolor', 'r');
h5 = plot(x(labels==5), RSS_cnmfe_norm(labels==5), 'o', 'markersize', 18, 'markerfacecolor', 'g');
temp = legend([h1, h2, h3, h4, h5], label_meanings, 'Interpreter', 'tex');
set(temp, 'fontsize', 20); 
xlabel('Operation #');
ylabel('normalized RSS');
set(gca, 'xtick', 2:2:10); 
set(gca, 'ytick', 1:0.01:1.02); 
ylim([0.999, 1.02]); 
xlim([0.5, 12]); 

if export_fig
    saveas(gcf, sprintf('%s/RSS.fig', output_folder));
    saveas(gcf, sprintf('%s/RSS.pdf', output_folder));
    saveas(gcf, sprintf('%s/RSS.png', output_folder));
end

%% save results
kt = 1;     % play one frame in every kt frames
save_avi = true;
t_begin = 1;
t_end = T;

center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2))*0.5;
range_res = [-1,1]*center_ac;
range_ac = center_ac+range_res+1;
multi_factor = 20;
center_Y = min(Y(:)) + multi_factor*center_ac;
range_Y = center_Y + range_res*multi_factor;
avi_filename =[video_folder, 'sim_snr1_decompose.avi'];
if ~exist(avi_filename, 'file')
    tmp_name = neuron.show_demixed_video(save_avi, kt, [t_begin, t_end], center_ac, range_ac, range_Y, multi_factor);
    movefile(tmp_name, avi_filename);
end
Ybg = neuron.reconstruct_background();
Ybg = neuron.reshape(Ybg, 1);
Y = neuron.reshape(Y,1);
Ysignal = Y - Ybg;

% show decomposed images for one example frame
figure('papersize', [d2+100, d1]/max(d1,d2)*5);
init_fig;
ind_frame = 578;
neuron.image(Y(:, ind_frame), [400, 2200]);
set(gca,'position',[0 0.03 .9 0.94],'units','normalized')
colorbar;
axis equal off tight;

hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end

figure('papersize', [d2+100, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ybg(:, ind_frame), [400, 2200]);
set(gca,'position',[0 0.03 .9 0.94],'units','normalized')
colorbar;
axis equal off tight;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_bg.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_bg.pdf', output_folder));
end

figure('papersize', [d2+100, d1]/max(d1,d2)*5);
init_fig;
ind_frame = 151;
neuron.image(neuron.A*neuron.C(:, ind_frame), [1, 81]);
colorbar;
set(gca,'position',[0 0.03 .9 0.94],'units','normalized')
axis equal off tight;
hold on;
if export_fig
    saveas(gcf, sprintf('%s/example_frame_ac.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_ac.pdf', output_folder));
end

figure('papersize', [d2+100, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ysignal(:, ind_frame) - neuron.A*neuron.C(:, ind_frame), [-40, 40]);
set(gca,'position',[0 0.03 .9 0.94],'units','normalized')
axis equal off tight; colorbar;
hold on;
if export_fig
    saveas(gcf, sprintf('%s/example_frame_res.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_res.pdf', output_folder));
end

figure('papersize', [d2+100, d1]/max(d1,d2)*5);
init_fig;
neuron.image(A*C(:, ind_frame), [1,81]); colorbar;
set(gca,'position',[0 0.03 .9 0.94],'units','normalized')
axis equal off tight;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_true_ac.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_true_ac.pdf', output_folder));
end

%% run CNMF and choose the one with the best performances
% initialization
Ain = A;
Cin = C;
res = neuron.reshape(Y, 1) - Ain*Cin;
res_ds = imresize(neuron.reshape(res,2), 0.5); 
[d1s, d2s, T] = size(res_ds); 
res_ds = reshape(res_ds, d1s*d2s, T); 
res_ds = imresize(res_ds, [size(res_ds,1), floor(T/2)]); 
tsub = 2;   % downsampling for fast initialization of BG
maxIter = 5;
neuron = neuron_basic; 
neuron.P.p = 2;
neuron.P.sn = neuron0.P.sn; 
try
    for m=1:16
        eval(sprintf('neuron_cnmf_%d=results_data.neuron_cnmf_%d;', m, m));
    end
catch
    for nb=1:16
        
        [bin, fin] = nnmf(res_ds, nb);
        %         fin = imresize(f, [nb, T]);
        bin = neuron.reshape(imresize(reshape(bin, d1s, d2s, []), [d1, d2]), 1);
        fin = imresize(fin, [size(fin,1), T]); 
%         bin = HALS_spatial(max(res, 0), bin, fin, [], maxIter);
%         fin = HALS_temporal(max(res, 0), bin, fin, maxIter);
        
        neuron_cnmf = neuron.copy();
        neuron_cnmf.A = Ain;
        neuron_cnmf.C = Cin;
        neuron_cnmf.b = bin;
        neuron_cnmf.f = fin;
        
        % iteratively update spatial and temporal components
        neuron_cnmf.P.mis_entries = [];
        for miter=1:2
            neuron_cnmf.updateSpatial(Y);
            neuron_cnmf.updateTemporal(Y);
        end
        neuron_cnmf.compress_results();
       
        eval(sprintf('neuron_cnmf_%d=neuron_cnmf.copy();', nb));
        clc;
        fprintf('number of background: %d\n', nb);
    end
    save(results_file, 'neuron_cnmf*', '-append');
end

%% evaluate performance of different rank
num_detected = zeros(1, 16);
sim_spatial = zeros(1, 16);
sim_temporal = zeros(1, 16);
figure('papersize', [10, 10]);
init_fig;

for nb=1:16
    eval(sprintf('neuron_cnmf = neuron_cnmf_%d.copy;', nb));
    A_cnmf = neuron_cnmf.A;
    C_cnmf = neuron_cnmf.C;
    match_cnmf = pair_neurons(A, C, A_cnmf, C_cnmf);
    ind_cnmf = ~isnan(match_cnmf.max_all);
    num_detected(nb) = sum(ind_cnmf);
    sim_spatial(nb) = mean(match_cnmf.max_spatial(ind_cnmf));
    sim_temporal(nb) = mean(match_cnmf.max_temporal(ind_cnmf));
    subplot(4,4,nb);
    
    plot(match_cnmf.max_spatial(ind_cnmf), match_cnmf.max_temporal(ind_cnmf), 'ok');
    box on;
    legend(sprintf('nb = %d', nb), 'location', 'northwest');
    %     xlabel('spatial similarity');
    %     ylabel('temporal similarity');
    set(gca, 'xtick', [0, 1]);
    set(gca, 'ytick', [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    if mod(nb, 4)==1
        ylabel('temporal');
    end
    if nb>12
        xlabel('spatial');
    end
end
% [~, nb_best] = max(sim_spatial+sim_temporal);
nb_best = 7; 
eval(sprintf('neuron_cnmf = neuron_cnmf_%d.copy;', nb_best));
save(results_file, 'neuron_cnmf', '-append');

if export_fig
    saveas(gcf, sprintf('%s/cnmf_nb_similarities.fig', output_folder));
    saveas(gcf, sprintf('%s/cnmf_nb_similarity.pdf', output_folder));
end


%% PCA/ICA
try
    A_ica_1 = results_data.A_ica_1;
    C_ica_1 = results_data.C_ica_1;
catch
    tic;
    [A_ica_1, C_ica_1] = run_pca_ica(neuron.reshape(Y,2), 250, 220);
    save(results_file, 'A_ica_1', 'C_ica_1', '-append');
    toc;
end

%% match neurons
% match CNMF-E initialization
match_cnmfe_init = pair_neurons(A, C, neuron_init.A, neuron_init.C_raw);
ind_cnmfe_init = ~isnan(match_cnmfe_init.max_all);

% match CNMF-E results
match_cnmfe = pair_neurons(A, C, neuron.A, neuron.C_raw);
ind_cnmfe = ~isnan(match_cnmfe.max_all);

% match PCA/ICA results
match_ica = pair_neurons(A, C, A_ica_1, C_ica_1);
ind_ica = ~isnan(match_ica.max_all);

% match CNMF results
A_cnmf = neuron_cnmf.A;
C_cnmf = neuron_cnmf.C;
match_cnmf = pair_neurons(A, C, A_cnmf, C_cnmf);
ind_cnmf = ~isnan(match_cnmf.max_all);


%% Figure 4C
figure('papersize', [6,6]);
init_fig;
set(gcf, 'defaultAxesfontsize', 25);
axes('position', [0.21, 0.17, 0.77, 0.8]);
hold on;
mksiz = 45;
scatter(match_ica.max_spatial(ind_ica), match_ica.max_temporal(ind_ica), 'sb','sizedata', mksiz);
scatter(match_cnmf.max_spatial(ind_cnmf), match_cnmf.max_temporal(ind_cnmf), 'vg',  'sizedata', mksiz);
scatter(match_cnmfe_init.max_spatial(ind_cnmfe_init), match_cnmfe_init.max_temporal(ind_cnmfe_init),...
    'dk',  'sizedata', mksiz);
scatter(match_cnmfe.max_spatial(ind_cnmfe), match_cnmfe.max_temporal(ind_cnmfe), ...
    'or',  'sizedata', mksiz);

box on;
% temp = legend(sprintf('PCA/ICA, n=%d', sum(ind_ica)), ...
%     sprintf('CNMF, n=%d', sum(ind_cnmf)),...
%     sprintf('CNMF-E initialization, n=%d', sum(ind_cnmfe_init)), ...
%     sprintf('CNMF-E final, n=%d', sum(ind_cnmfe)), 'location', 'southwest');
temp = legend('PCA/ICA',  'CNMF', 'CNMF-E initialization','CNMF-E', 'location', 'southwest');

set(temp, 'fontsize', 18);
xlabel('Spatial similarity');
ylabel('Temporal similarity');
xlim([0.9, 1]);
ylim([0.402, 1]);
set(gca, 'xtick', 0.9:0.05:1);
set(gca, 'ytick', 0.4:0.2:1);
apply_fisher = true;

if export_fig
    saveas(gcf, sprintf('%s/example_sim_similarity_with_cnmf.fig', output_folder));
    saveas(gcf, sprintf('%s/example_sim_similarity_with_cnmf.pdf', output_folder));
end

if apply_fisher
    temp = get(gca, 'children');
    max_x = 0.0;
    max_y = 0.0;
    for i=1:4
        xdata = get(temp(i), 'Xdata');
        ydata = get(temp(i), 'Ydata');
        max_x = max(max_x, max(xdata));
        max_y = max(max_y, max(ydata));
        set(temp(i), 'Xdata', atanh(xdata));
        set(temp(i), 'Ydata', atanh(ydata));
    end
    xlim(atanh([-0.5, max_x]));
    ylim(atanh([-0.6, max_y]));
    v_xtick = [0, 0.6, 0.9, 0.99, 0.999];
    v_ytick = [0, 0.5, 0.8, 0.95, 0.99];
    set(gca, 'xtick', atanh(v_xtick));
    set(gca, 'xticklabel', v_xtick);
    set(gca, 'ytick', atanh(v_ytick));
    set(gca, 'yticklabel', v_ytick);
    
    if export_fig
        saveas(gcf, sprintf('%s/example_sim_similarity_fisher.fig', output_folder));
        saveas(gcf, sprintf('%s/example_sim_similarity_fisher.pdf', output_folder));
    end
    
end

%% correlation structure. Figure 4D
close all;
C_true = corr(C') - eye(size(C, 1));
figure('papersize', [9, 9]);
init_fig;
set(gcf, 'defaultaxesfontsize', 30);
v_alpha = 1;


% pca/ica, thresholded.
axes('position', [0.57, 0.57, 0.41, 0.41]);
axis([-0.2, 1, -0.2, 1]); hold on;
set(gca, 'xtick', 0:0.5:1);
set(gca, 'ytick', 0:0.5:1);
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
box on;
temp_true = C_true(ind_ica, ind_ica);
tmp_C_ica = C_ica_1;
tmp_C_ica(bsxfun(@lt, C_ica_1, 3*std(C_ica_1, [], 2))) = 0;
C_pca_ica = corr(tmp_C_ica(match_ica.ind_max(ind_ica), :)') - eye(sum(ind_ica));
% plot(temp_true(:), C_pca_ica(:), '*b', 'markersize', 5);  hold on;
% scatter(temp_true(:), C_pca_ica(:), 'v', 'filled','markerfacecolor', 'm', ...
scatter(temp_true(:), C_pca_ica(:), 'pm');
plot([-0.2, 1], [-0.2, 1], 'k');
text(-0.15, 0.93, 'PCA/ICA', 'fontsize', 30, 'color', 'k')
text(-0.15, 0.8, '+thresholding', 'fontsize', 30, 'color', 'k')
% temp = legend('PCA/ICA+thresholding', 'location', 'northwest');
% set(temp, 'fontsize', 20);

% pca/ica, without thresholding
C_pca_ica = corr(C_ica_1(match_ica.ind_max(ind_ica), :)') - eye(sum(ind_ica));
axes('position', [0.15, 0.57, 0.41, 0.41]);
axis([-0.2, 1, -0.2, 1]); hold on;
set(gca, 'xtick', 0:0.5:1);
set(gca, 'ytick', 0:0.5:1);
set(gca, 'xticklabel', []);
box on;
ylabel('Estimated corr.');
% scatter(temp_true(:), C_pca_ica(:), 's', 'filled','markerfacecolor', 'b', ...
%     'markerfacealpha', v_alpha);
scatter(temp_true(:), C_pca_ica(:), 'sb');
plot([-0.2, 1], [-0.2, 1], 'k');
text(-0.15, 0.93, 'PCA/ICA', 'fontsize', 30, 'color', 'k')
% temp = legend('PCA/ICA', 'location', 'northwest');
% set(temp, 'fontsize', 20);

% cnmf
temp_true = C_true(ind_cnmf, ind_cnmf);
C_cnmf = corr(neuron_cnmf.C(match_cnmf.ind_max(ind_cnmf), :)') - eye(sum(ind_cnmf));
axes('position', [0.15, 0.15, 0.41, 0.41]);
axis([-0.2, 1, -0.2, 1]); hold on;
set(gca, 'xtick', 0:0.5:1);
set(gca, 'ytick', 0:0.5:1);
box on;
% scatter(temp_true(:), C_cnmf(:), 'v', 'filled','markerfacecolor', 'g', ...
%     'markerfacealpha', v_alpha);
scatter(temp_true(:), C_cnmf(:), 'vg');
diff_corr_cnmf = C_cnmf(:)-temp_true(:);
ylabel('Estimated corr.');
xlabel('True corr.');
plot([-0.2, 1], [-0.2, 1], 'k');
text(-0.15, 0.93, 'CNMF', 'fontsize', 30, 'color', 'k')
% temp =legend('CNMF', 'location', 'northwest');
% set(temp, 'fontsize', 20);

% cnmfe
temp_true = C_true(ind_cnmfe, ind_cnmfe);
C_cnmfe = corr(neuron.C_raw(match_cnmfe.ind_max(ind_cnmfe), :)') - eye(sum(ind_cnmfe));
axes('position', [0.57, 0.15, 0.41, 0.41]);
axis([-0.2, 1, -0.2, 1]); hold on;
set(gca, 'xtick', 0:0.5:1);
set(gca, 'ytick', 0:0.5:1);
set(gca, 'yticklabel', []);
% box on; scatter(temp_true(:), C_cnmfe(:), 'o', 'filled','markerfacecolor', 'r', ...
%     'markerfacealpha', v_alpha);
box on; scatter(temp_true(:), C_cnmfe(:), 'or');
diff_corr_cnmfe = C_cnmfe(:) - temp_true(:);
axis tight;
plot([-0.2, 1], [-0.2,1], '-k');
xlim([-0.2, 1]);
ylim([-0.2, 1]);
box on;
xlabel('True corr.');
plot([-0.2, 1], [-0.2, 1], 'k');
text(-0.15, 0.93, 'CNMF-E', 'fontsize', 30, 'color', 'k')
% temp = legend('CNMF-E', 'location', 'northwest');
% set(temp, 'fontsize', 20);
if export_fig
    saveas(gcf, sprintf('%s/example_sim_correlation_separate.fig', output_folder));
    %         saveas(gcf, sprintf('%s/example_sim_correlation_separate.pdf', output_folder));
    saveas(gcf, sprintf('%s/example_sim_correlation_separate.png', output_folder));
end

%% distribution of the error in recovering the true corr. coeff.
figure('papersize', [6, 6]);
init_fig;
bins = linspace(-0.5, 0.5, 100);
c1 = hist(diff_corr_cnmf, bins);
c2 = hist(diff_corr_cnmfe, bins);
fill([bins, bins(1)], [c1, c1(1)]/sum(c1), 'g', 'FaceAlpha', 0.3, 'edgecolor', 'g');
hold on
fill([bins, bins(1)], [c2, c2(1)]/sum(c2), 'r', 'FaceAlpha', 0.3, 'edgecolor', 'r')
xlim([-0.3, 0.3]);
xlabel('\Delta corr. coeff.');
ylabel('Frequency');
legend('CNMF', 'CNMF-E');
if export_fig
    saveas(gcf, sprintf('%s/sim_delta_corr.fig', output_folder));
    %     saveas(gcf, sprintf('%s/example_sim_correlation_separate.pdf', output_folder));
    saveas(gcf, sprintf('%s/sim_delta_corr.pdf', output_folder));
end

%% select overlapped neurons as examples. Figure 4A
% cell_ids = [55, 106, 143];
cell_ids = [1, 151,187];
cell_ids_cnmfe = match_cnmfe.ind_max(cell_ids);
cell_ids_ica = match_ica.ind_max(cell_ids);
cell_ids_cnmf = match_cnmf.ind_max(cell_ids);
% cell_ids_cnmfe(isnan(cell_ids_cnmfe)) = match_cnmfe.ind_spatial(cell_ids(isnan(cell_ids_cnmfe) ));
% cell_ids_ica(isnan(cell_ids_ica)) = match_ica.ind_spatial(cell_ids(isnan(cell_ids_ica)));

temp = neuron.reshape(sum(A(:, cell_ids), 2), 2);
[tmp_r, tmp_c] = find(temp>1e-5);
xmin = min(tmp_c);
xmax = max(tmp_c);
ymin = min(tmp_r);
ymax = max(tmp_r);
x_center = round(xmin/2+xmax/2);
y_center = round(ymin/2+ymax/2);
ctr = sub2ind([d1,d2], y_center, x_center);

% spatial components
figure('papersize', [(xmax-xmin+1),(ymax-ymin+1)]*0.1);
init_fig;
img = neuron.reshape(A(:, cell_ids), 2);
axes('position', [0, 0, 1, 1]);
imagesc(img/max(img(:)));
axis equal off;
axis([xmin, xmax, ymin, ymax]);
if export_fig
    saveas(gcf, sprintf('%s/example_spatial_true.fig', output_folder));
    saveas(gcf, sprintf('%s/example_spatial_true.pdf', output_folder));
end

% temporal components
figure('papersize', [12, 2]);
init_fig;
axes('position', [0.05, 0.05, 0.9, 0.9]); hold on;
y = C(cell_ids(1), :);
plot(y/max(y)+1.6, 'r', 'linewidth', 3);

y = C(cell_ids(2), :);
plot(y/max(y)+0.8, 'g', 'linewidth', 3);

y = C(cell_ids(3), :);
plot(y/max(y) , 'b', 'linewidth', 3);
axis tight;
axis off;
xlim([0, 1200]);
if export_fig
    saveas(gcf, sprintf('%s/example_temporal_true.fig', output_folder));
    saveas(gcf, sprintf('%s/example_temporal_true.pdf', output_folder));
end

% pca/ica
A_ica = A_ica_1;
C_ica = C_ica_1;
figure('papersize', [(xmax-xmin+1),(ymax-ymin+1)]*0.1);
init_fig;
img = neuron.reshape(A_ica(:, max(1,cell_ids_ica)), 2);
for m=1:3
    if isnan(cell_ids_ica(m))
        img(:, :, m) = 0;
    end
end
img(img<0) = 0;
imagesc(img/max(img(:)));
axis equal off;
set(gca, 'position', [0, 0, 1, 1]); hold on;

axis([xmin, xmax, ymin, ymax]);
if export_fig
    saveas(gcf, sprintf('%s/example_spatial_ica.fig', output_folder));
    saveas(gcf, sprintf('%s/example_spatial_ica.pdf', output_folder));
    
    axis tight;
    saveas(gcf, sprintf('%s/example_spatial_ica_whole.fig', output_folder));
    saveas(gcf, sprintf('%s/example_spatial_ica_whole.pdf', output_folder));
    
end

figure('papersize', [12, 2]);
init_fig;
axes('position', [0.05, 0.05, 0.9, 0.9]); hold on;
if ~isnan(cell_ids_ica(1))
    y = C_ica(cell_ids_ica(1), :);
    plot(y/max(y)+1.6, 'r', 'linewidth', 3);
end

if ~isnan(cell_ids_ica(2))
    y = C_ica(cell_ids_ica(2), :);
    plot(y/max(y)+0.8, 'g', 'linewidth', 3);
end

if ~isnan(cell_ids_ica(3))
    y = C_ica(cell_ids_ica(3), :);
    plot(y/max(y) , 'b', 'linewidth', 3);
end
axis tight;
xlim([0, 1200]);
axis off;
if export_fig
    saveas(gcf, sprintf('%s/example_temporal_ica.fig', output_folder));
    saveas(gcf, sprintf('%s/example_temporal_ica.pdf', output_folder));
end

% CNMF-E
figure('papersize', [(xmax-xmin+1),(ymax-ymin+1)]*0.1);
init_fig;
axes('position', [0, 0, 1, 1]);
img = neuron.reshape(neuron.A(:, cell_ids_cnmfe), 2);
imagesc(img/max(img(:)));

axis equal off;
axis([xmin, xmax, ymin, ymax]);
if export_fig
    saveas(gcf, sprintf('%s/example_spatial_cnmfe.fig', output_folder));
    saveas(gcf, sprintf('%s/example_spatial_cnmfe.pdf', output_folder));
end
%
figure('papersize', [12,2]);
init_fig;
axes('position', [0.05, 0.05, 0.9, 0.9]); hold on;
% y1 = neuron.C_raw(cell_ids_cnmfe(1), :);
y = neuron.C(cell_ids_cnmfe(1), :);
% plot(y1/max(y1)+1.6, 'k', 'linewidth', 3);
plot(y/max(y)+1.6, 'r', 'linewidth', 3);

% y1 = neuron.C_raw(cell_ids_cnmfe(2), :);
y = neuron.C(cell_ids_cnmfe(2), :);
% plot(y1/max(y1)+0.8, 'k', 'linewidth', 3);
plot(y/max(y)+0.8, 'g', 'linewidth', 3);

% y1 = neuron.C_raw(cell_ids_cnmfe(3), :);
y = neuron.C(cell_ids_cnmfe(3), :);
% plot(y1/max(y1), 'k', 'linewidth', 3);
plot(y/max(y) , 'b', 'linewidth', 3);
axis tight;
xlim([0, 1200]);
axis off;
if export_fig
    saveas(gcf, sprintf('%s/example_temporal_cnmfe.fig', output_folder));
    saveas(gcf, sprintf('%s/example_temporal_cnmfe.pdf', output_folder));
end


% CNMF
figure('papersize', [(xmax-xmin+1),(ymax-ymin+1)]*0.1);
init_fig;
axes('position', [0, 0, 1, 1]);
img = neuron_cnmf.reshape(neuron_cnmf.A(:, cell_ids_cnmf), 2);
imagesc(img/max(img(:)));

axis equal off;
axis([xmin, xmax, ymin, ymax]);
if export_fig
    saveas(gcf, sprintf('%s/example_spatial_cnmf.fig', output_folder));
    saveas(gcf, sprintf('%s/example_spatial_cnmf.pdf', output_folder));
end
%
figure('papersize', [12,2]);
init_fig;
axes('position', [0.05, 0.05, 0.9, 0.9]); hold on;
% y1 = neuron.C_raw(cell_ids_cnmfe(1), :);
y = neuron_cnmf.C(cell_ids_cnmf(1), :);
% plot(y1/max(y1)+1.6, 'k', 'linewidth', 3);
plot(y/max(y)+1.6, 'r', 'linewidth', 3);

% y1 = neuron.C_raw(cell_ids_cnmfe(2), :);
y = neuron_cnmf.C(cell_ids_cnmf(2), :);
% plot(y1/max(y1)+0.8, 'k', 'linewidth', 3);
plot(y/max(y)+0.8, 'g', 'linewidth', 3);

% y1 = neuron.C_raw(cell_ids_cnmfe(3), :);
y = neuron_cnmf.C(cell_ids_cnmf(3), :);
% plot(y1/max(y1), 'k', 'linewidth', 3);
plot(y/max(y) , 'b', 'linewidth', 3);
axis tight;
xlim([0, 1200]);
axis off;
if export_fig
    saveas(gcf, sprintf('%s/example_temporal_cnmf.fig', output_folder));
    saveas(gcf, sprintf('%s/example_temporal_cnmf.pdf', output_folder));
end

%% check the dependence on noise levels
noise_amp = 1:6;
min_corr_vec = [0.9, 0.8, 0.8, 0.8, 0.6, 0.6];
min_pnr_vec = [15, 10, 10, 8, 6, 6];
A_ica = zeros(d1*d2, 300, 6);
C_ica = zeros(300, T, 6);
nPCs_vec = [600, 600, 600, 600, 600, 600];
nICs_vec = [300, 300, 300, 300, 300, 300];
save temp noise_amp;
try
    for m=1:6
        eval(sprintf('%s=results_data.%s; ', nam));
        eval(sprintf('%s_init=results.%s_init; ', nam));
    end
    A_ica = results_data.A_ica;
    C_ica = results_data.C_ica;
catch
    %%
    for m_noise =1:6
        nam = sprintf('neuron_%d', noise_amp(m_noise));
        if exist(nam, 'var')
            continue;
        end
        disp(m_noise);
        %% simulate data
        Y = A*C + bsxfun(@plus, Bf, Bc) + noise_amp(m_noise)*E;  % observed fluorescence
        Y = round(Y);
        %% run PCA/ICA method
        [A_ica(:, :, m_noise), C_ica(:, :, m_noise)] = run_pca_ica(neuron.reshape(Y,2),...
            nPCs_vec(m_noise), nICs_vec(m_noise));
        
        %% run CNMF-E
        Y = neuron.reshape(round(Y), 2);
        Ysiz = size(Y);
        nam_data = sprintf('%sdata_snr_%d.mat', results_folder, m_noise);
        if ~exist(nam_data, 'file')
            save(nam_data, 'Y', 'Ysiz', '-v7.3');
        end
        neuron.select_data(nam_data);
        Y = neuron.reshape(Y,1);
        rng(1);
        neuron.options.min_corr = min_corr_vec(m_noise);
        neuron.options.min_pnr = min_pnr_vec(m_noise);
        if m_noise>3
            neuron.options.deconv_options.smin = -2;
        end
        neuron.getReady(pars_envs);
        [center, Cn, PNR] = neuron.initComponents_parallel(K, [], save_initialization, use_parallel);
        neuron.merge_neurons_dist_corr();
        
        %% estimate the background components
        neuron.update_background_parallel(use_parallel);
        neuron_init = neuron.copy();
        
        %% pick neurons from the residual
        neuron.initComponents_residual_parallel([], save_initialization, use_parallel);
        neuron.merge_neurons_dist_corr();
        
        %% correct the inaccurate extraction of neural traces
        use_c_hat = false;
        neuron.update_temporal_parallel(use_parallel, use_c_hat);
        neuron.update_spatial_parallel(use_parallel);
        neuron.merge_neurons_dist_corr();
        neuron.remove_false_positives(); 
        
        %% udpate spatial&temporal components,
        for m=1:3
            if mod(m, 2) ==0
                neuron.update_background_parallel(use_parallel);
            else
                % update temporal
                neuron.update_temporal_parallel(use_parallel);
                % update _spatial
                neuron.update_spatial_parallel(use_parallel);
            end
        end
        
        eval(sprintf('%s=neuron.copy; ', nam));
        eval(sprintf('%s_init=neuron_init.copy; ', nam));
        eval(sprintf('save temp.mat %s_init %s -append', nam, nam));
        
        %% 
        temp = pair_neurons(A, C, neuron.A, neuron.C); 
        if sum(isnan(temp.ind_max))~=0
            disp(m_noise); 
            pause; 
        end
        
        %%
        if m_noise==6
            avi_filename = [video_folder, 'sim_snr6_decompose.avi'];
            if ~exist(avi_filename, 'file')
                tmp_name = neuron.show_demixed_video(save_avi, kt, [t_begin, t_end], center_ac, range_ac, range_Y, multi_factor);
                movefile(tmp_name, avi_filename);
            end
        end
    end
    delete temp.mat;
%     save(results_file, 'neuron_*', 'A_ica', 'C_ica', '-append');
end

%% calculate the quality of the extracted neurons
temporal_init_raw = nan(200, 6);
spatial = nan(200, 6);
temporal = nan(200, 6);
temporal_raw = nan(200,6);
spatial_ica = nan(200,6);
temporal_ica = nan(200, 6);

C_norm = bsxfun(@times, C, 1./sqrt(sum(C.^2, 2)));

for m=1:6
    eval(sprintf('neuron=neuron_%d;', m));
    eval(sprintf('neuron_init=neuron_%d_init;', m));
    %     eval(sprintf('neuron_cnmf=neuron_cnmf_noise_%d;', m));
    
    % match PCA/ICA results
    match_ica = pair_neurons(A, C, squeeze(A_ica(:, :, m)), squeeze(C_ica(:, :, m)));
    ind_ica = (~isnan(match_ica.ind_max));
    spatial_ica(ind_ica, m) = match_ica.max_spatial(ind_ica);
    temporal_ica(ind_ica, m) = match_ica.max_temporal(ind_ica);
    
    % match CNMF-E results
    match_cnmfe = pair_neurons(A, C, neuron.A, neuron.C);
    ind_cnmfe = (~isnan(match_cnmfe.ind_max));
    spatial(ind_cnmfe,m) = match_cnmfe.max_spatial(ind_cnmfe);
    temporal(ind_cnmfe,m) = match_cnmfe.max_temporal(ind_cnmfe);
    temp = neuron.C_raw(match_cnmfe.ind_spatial, :);
    temp = bsxfun(@times, temp, 1./sqrt(sum(temp.^2, 2)));
    temp2 = sum(temp.*C_norm, 2);
    temporal_raw(ind_cnmfe, m) = temp2(ind_cnmfe);
    
    %     % match CNMF results
    %     match_cnmf = pair_neurons(A, C, neuron_cnmf.A, neuron_cnmf.C);
    %     ind_cnmf = (~isnan(match_cnmf.ind_max));
    %     spatial_cnmf(ind_cnmf,m) = match_cnmf.max_spatial(ind_cnmf);
    %     temporal_cnmf(ind_cnmf,m) = match_cnmf.max_temporal(ind_cnmf);
end

save similarity_comparison spatial_* temporal_* spatial temporal;

%% number of false negatives. Figure 4D
figure('papersize', [5.8,6]);
init_fig;
set(gcf, 'defaultaxesfontsize', 30);
axes('position', [0.3, 0.22, 0.6875, 0.77]);
hold on;
plot(noise_amp(1:6), sum(isnan(spatial_ica), 1), '-sb', 'markersiz', 15);
hold on;
% plot(noise_amp(1:6), sum(isnan(spatial_init), 1), '-*g', 'markersize', 10);
plot(noise_amp(1:6), sum(isnan(spatial), 1), '-or', 'markersize', 15);
axis tight;
ylim([-10, 180]);
set(gca, 'xtick', 1:6);
set(gca, 'ytick', 0:50:200);
xlabel('SNR reduction factor');
ylabel('# of missed neurons');
box on;
temp = legend('PCA/ICA', 'CNMF-E', 'location', 'northwest');
set(temp, 'fontsize', 23); 
xlim([0.5, 6.5]);

if export_fig
    saveas(gcf, sprintf('%s/false_negatives.fig', output_folder));
    saveas(gcf, sprintf('%s/false_negatives.pdf', output_folder));
end

%% show spatial similarities. Figure 4E
figure('papersize', [5.8,6]);
init_fig;
set(gcf, 'defaultaxesfontsize', 30);
axes('position', [0.3, 0.22, 0.6875, 0.77]);
hold on;

mean_ica = nanmean(spatial_ica);
std_ica = nanstd(spatial_ica);
mean_cnmfe = nanmean(spatial);
std_cnmfe = nanstd(spatial);

plot(1:6, mean_ica, '-sb', 'markersize', 15); hold on;
% plot(1:6, mean_cnmf, '-*g', 'markersize', 15);
% errorbar(1:6, mean_init,std_init, '-*g');
plot(1:6, mean_cnmfe, '-or', 'markersize', 15);
set(gca, 'xtick', 1:6);
set(gca, 'ytick', 0.5:0.05:1);
ylim([0.89, 1.001]);
xlim([0.5, 6.5]);
xlabel('SNR reduction factor');
ylabel('Spatial similarity');
box on;
temp = legend('PCA/ICA', 'CNMF-E', 'location', 'southwest');
set(temp, 'fontsize', 21); 
if export_fig
    saveas(gcf, sprintf('%s/noise_spatial_similarity.fig', output_folder));
    saveas(gcf, sprintf('%s/noise_spatial_similarity.pdf', output_folder));
end

%% show temporal similarities. Figure 4F
figure('papersize', [5.8,6]);
init_fig;
set(gcf, 'defaultaxesfontsize', 30);
axes('position', [0.3, 0.22, 0.6875, 0.77]);
hold on;

mean_ica = nanmean(temporal_ica);
std_ica = nanstd(temporal_ica);
mean_raw = nanmean(temporal_raw);
std_raw = nanstd(temporal_raw);
mean_cnmfe = nanmean(temporal);
std_cnmfe = nanstd(temporal);

plot(1:6, mean_ica, '-sb', 'markersize', 15);
% plot(1:6, mean_cnmf, '-*m', 'markersize', 15);
% errorbar(1:6, mean_init,std_init, '-*g');
plot(1:6, mean_raw, '-<c', 'markersize', 15);
plot(1:6, mean_cnmfe, '-or', 'markersize', 15);
set(gca, 'xtick', 1:6);
set(gca, 'xticklabel', 1:6);
set(gca, 'ytick', 0.2:0.2:1);
ylim([0.35, 1.01]);
xlim([0.5, 6.5]);

xlabel('SNR reduction factor');
ylabel('Temporal similarity');
box on;
temp = legend('PCA/ICA', 'CNMF-E, no denoising', 'CNMF-E', 'location', 'southwest');
set(temp, 'fontsize', 21); 
if export_fig
    saveas(gcf, sprintf('%s/noise_temporal_similarity.fig', output_folder));
    saveas(gcf, sprintf('%s/noise_temporal_similarity.pdf', output_folder));
end

%% use one example to show the extracted traces. Figure 4G
ind = find(~isnan(spatial_ica(:, 6)));
n = 2;
ind_ica = match_ica.ind_max(ind(n));
ind_cnmfe = match_cnmfe.ind_max(ind(n));

figure('papersize', [5.5, 4]);
init_fig;
y_true = C(ind(n), :);
y_ica = squeeze(C_ica(ind_ica, :, 6));
y_raw = neuron_6.C_raw(ind_cnmfe, :);
y_cnmfe = neuron_6.C(ind_cnmfe, :);

axes('position', [0.01, 0, 0.98, 1]);
hold on;
plot(y_true/max(y_true)+2.5, 'k');

plot(y_ica/max(y_ica)+1.25, 'b');
text(300, 2.25, sprintf('corr=%.3f', match_ica.max_temporal(ind(n))), ...
    'color', 'b', 'fontweight', 'bold', 'fontsize', 28);

plot(y_raw/max(y_raw), 'c');

plot(y_cnmfe/max(y_cnmfe), 'r');
text(300, 0.76, sprintf('corr=%.3f', temporal_raw(ind(n), 6)), ...
    'color', 'c', 'fontweight', 'bold', 'fontsize', 28);
text(300, 0.5, sprintf('corr=%.3f', temporal(ind(n), 6)), ...
    'color', 'r', 'fontweight', 'bold', 'fontsize', 28);
axis off; axis tight;
ylim([-0.25, 3.6]);
xlim([0, 1200]); 
% legend('true', 'PCA/ICA', 'CNMF-E, no denoising', 'CNMFE', 'orientation',...
%     'horizental');
set(gca, 'fontsize', 30);

if export_fig
    saveas(gcf, sprintf('%s/noise_traces.fig', output_folder));
    saveas(gcf, sprintf('%s/noise_traces.pdf', output_folder));
end

%% save results



