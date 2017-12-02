%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_HIPPOCAMPUS_subfigs'];

results_folder = sprintf('%sfig_HIPPOCAMPUS%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_HIPPOCAMPUS_results.mat'];
cellsort_folder = [code_folder, 'CellSort'];
addpath(cellsort_folder);
demixed_video = [video_folder, filesep, 'HIPPOCAMPUS_decompose.avi'];
intervention_video = [video_folder, filesep, 'intervention.avi'];
% overlap_video = [video_folder, filesep, 'PFC_overlapping.avi'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file, 'Writable', true);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% run CNMF-E or load results
nam = get_fullname('../../../data/bma22_epm_motioncorrected_round1_cropped_correction.mat');
try
    cnmfe_path = results_data.cnmfe_path;
catch
    cnmfe_path = [];
end
if exist(cnmfe_path, 'file')
    load(cnmfe_path);
else
    run_cnmfe;
    results_data.cnmfe_path = cnmfe_path;
end
frame_range = neuron.frame_range;
if ~exist('Y', 'var')
    Y = neuron.load_patch_data([], frame_range);
end

% load raw data and reconstruct background
Y = neuron.reshape(Y, 1);
Ybg = neuron.reconstruct_background();
Ybg = neuron.reshape(Ybg, 1);
Ysignal = neuron.reshape(double(Y), 1) - neuron.reshape(Ybg, 1);

%% run PCA/ICA analysis or load results
nPCs = 100;
nICs = 50;
try
    neuron_ica = results_data.neuron_ica;
    A_ica_before_trim = results_data.A_ica_before_trim;
catch
    run_ica;
    
    %%
    results_data.A_ica_before_trim = A_ica_before_trim;
    results_data.A_ica = A_ica;
    results_data.C_ica = C_ica;
    results_data.time_cost = time_cost;
    
    neuron_ica.compress_results();
    results_data.neuron_ica = neuron_ica;
end

%% compute SNRs
K_cnmfe = size(neuron.C_raw, 1);
K_ica = size(neuron_ica.C, 1);

snr_ica = var(neuron_ica.C, 0, 2)./var(neuron_ica.C_raw-neuron_ica.C, 0, 2);
snr_cnmfe = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);

[~, srt] = sort(snr_cnmfe, 'descend');
neuron.orderROIs(srt);
snr_cnmfe = snr_cnmfe(srt);

[~, srt] = sort(snr_ica, 'descend');
neuron_ica.orderROIs(srt);
A_ica_before_trim = A_ica_before_trim(:, srt);
snr_ica = snr_ica(srt);

%% match neurons
match_cnmfe = pair_neurons(neuron.A, neuron.C, neuron_ica.A, neuron_ica.C_raw);

ids_match = find(~isnan(match_cnmfe.ind_max));
ids_ica = match_cnmfe.ind_max(ids_match);
K_match = length(ids_match);

ind = 1:size(neuron.A,2);
ind(ids_match) = [];
srt = [ids_match, ind];
neuron.orderROIs(srt);
snr_cnmfe = snr_cnmfe(srt);

ind = 1:size(neuron_ica.A, 2);
ind(ids_ica) = [];
srt = [ids_ica, ind];
neuron_ica.orderROIs(srt);
snr_ica = snr_ica(srt);
A_ica_before_trim = A_ica_before_trim(:, srt);

%% match neurons
match_cnmfe = pair_neurons(neuron.A, neuron.C, neuron_ica.A, neuron_ica.C_raw);
match_ica = pair_neurons(neuron_ica.A, neuron_ica.C_raw, neuron.A, neuron.C_raw);
snr_ica = var(neuron_ica.C, 0, 2)./var(neuron_ica.C_raw-neuron_ica.C, 0, 2);
snr_cnmfe = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);

%% compute correlation images
d1 = neuron.options.d1;
d2 = neuron.options.d2;
create_correlation_images;

%%  contours of the detected neurons, CNMF-E
neuron.Cn = Cn_nobg;
neuron.PNR = PNR_filter;
if ~exist('Coor_cnmfe', 'var') || (length(Coor_cnmfe)~=size(neuron.A, 2))
    Coor_cnmfe = neuron.get_contours(0.8);
    neuron.Coor = Coor_cnmfe;
end
if ~exist('Coor_ica', 'var')|| (length(Coor_ica)~=size(neuron_ica.A, 2))
    Coor_ica = neuron_ica.get_contours(0.8);
    neuron_ica.Coor = Coor_ica;
end

neuron_cnmfe = neuron.copy(); 
%%

%% plot contours
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;

imagesc(Cn_nobg, [0, 1]); colormap gray; hold on; %colorbar;
axis equal off tight;
set(gca,'position',[0.01 0.03 .85 0.9],'units','normalized');
for m=1:length(Coor_cnmfe)
    temp = Coor_cnmfe{m};
    if m==1
        h1 = fill(temp(1,:), temp(2,:), 'r', 'facealpha', 0.3,'edgecolor', 'w');
    else
        fill(temp(1,:), temp(2,:), 'r', 'facealpha', 0.3,'edgecolor', 'w');
    end
    [~, ind] = max(temp(1,:)); %randi(length(temp), 1);
    text(temp(1,ind), temp(2,ind), num2str(m), 'color', 'r', ...
        'fontsize', 14,'fontweight', 'bold');
end

for m=1:length(Coor_ica)
    temp = medfilt1(Coor_ica{m}')';
    temp = temp(:, 3:end);
    if m==1
        h2 = fill(temp(1,:), temp(2,:), 'g', 'facealpha', 0.3,'edgecolor', 'w');
    else
        fill(temp(1,:), temp(2,:), 'g', 'facealpha', 0.3,'edgecolor', 'w');
    end
    [~, ind] = min(temp(1,:)); %randi(length(temp), 1);
    text(temp(1,ind)-7, temp(2,ind), num2str(m), 'color', 'g', ...
        'fontsize', 14, 'fontweight', 'bold');
end

colorbar; 
set(gca, 'fontsize', 14);
legend([h1, h2], {'CNMF-E', 'PCA/ICA'}, 'location', 'southeast');
if export_fig
    saveas(gcf, sprintf('%s/contours.fig', output_folder));
    saveas(gcf, sprintf('%s/contours.pdf', output_folder));
end

%% spatial components of all matched neurons
figure('papersize', [7, 4]);
init_fig;

ctr = neuron.estCenter();
dw = 1/7;
dh = 1/4;
w = 3*neuron.options.gSiz + 1;
K_cnmfe = size(neuron.A,2);
K_ica = size(neuron_ica.A, 2);
K_match = length(ids_match);

% matched neuron
for m=1:K_cnmfe
    if m<=K_match
        ii = mod(m, 7);
        jj = ceil(m/7);
    else
        ii = m-K_match;
        jj = 4;
    end
    axes('position', [(ii-1)*dw+1*(ii==0)+0.001,(4-jj)*dh+0.001,dw-0.002, dh-0.001]);
    
    r0 = round(ctr(m,1));
    c0 = round(ctr(m,2));
    if (r0-w)<1
        r0 = w+1;
    elseif (r0+w)>d1
        r0 = d1-w;
    end
    
    if (c0-w)<1
        c0 = w+1;
    elseif (c0+w)>d2
        c0 = d2-w;
    end
    
    img = neuron.reshape(neuron.A(:, m), 2);
    temp = img(r0+(-w:w), c0+(-w:w));
    imagesc(temp/max(temp(:)), [0,1]);
    axis equal off; hold on;
    text(1, 15, num2str(m), 'color', 'w','fontweight', 'bold', 'fontsize', 20);
    if m==K_match+1
        plot([2*w+1,1,1,2*w+1], [1, 1,2*w+1,2*w+1], 'r', 'linewidth', 2);
    elseif m>K_match && m<K_cnmfe
        plot([2*w+1,-5], [2*w+1,2*w+1], 'r', 'linewidth', 2);
        plot([2*w+1,-5], [1,1], 'r', 'linewidth', 2);
    elseif m==K_cnmfe
        plot([-5, 2*w+1, 2*w+1, -5], [2*w+1, 2*w+1, 1, 1], 'r', 'linewidth', 2);
    end
    xlim([0, 2*w+2]); ylim([0, 2*w+2]);
end


if export_fig
    saveas(gcf, sprintf('%s/cnmfe_spatial.fig', output_folder) );
    saveas(gcf, sprintf('%s/cnmfe_spatial.pdf', output_folder));
end

%% show ICA neurons
figure('papersize', [7, 4]);
init_fig;

ctr = neuron_ica.estCenter();

% matched neuron
for m=1:K_ica
    if m<=K_match
        ii = mod(m, 7);
        jj = ceil(m/7);
    else
        ii = m-K_match;
        jj = 4;
    end
    axes('position', [(ii-1)*dw+1*(ii==0)+0.001,(4-jj)*dh+0.001,dw-0.002, dh-0.001]);
    
    r0 = round(ctr(m,1));
    c0 = round(ctr(m,2));
    if (r0-w)<1
        r0 = w+1;
    elseif (r0+w)>d1
        r0 = d1-w;
    end
    
    if (c0-w)<1
        c0 = w+1;
    elseif (c0+w)>d2
        c0 = d2-w;
    end
    
    img = neuron.reshape(A_ica_before_trim(:, m), 2);
    temp = img(r0+(-w:w), c0+(-w:w));
    imagesc(temp/max(temp(:)), [-0.3, 1]);
    axis equal off; hold on;
    text(1, 15, num2str(m), 'color', 'w','fontweight', 'bold', 'fontsize', 20);
    xlim([0, 2*w+1]); ylim([0, 2*w+1]);
    if jj<4
        continue;
    end
    if m==K_match+1
        plot([2*w+1,0,0,2*w+1], [0, 0,2*w+1,2*w+1], 'r', 'linewidth', 2);
    elseif m>K_match && m<K_ica
        plot([2*w+1,-5], [2*w+1,2*w+1], 'r', 'linewidth', 2);
        plot([2*w+1,-5], [0,0], 'r', 'linewidth', 2);
    elseif m==K_ica
        plot([-5, 2*w+1, 2*w+1, -5], [2*w+1, 2*w+1, 0, 0], 'r', 'linewidth', 2);
    end
end


if export_fig
    saveas(gcf, sprintf('%s/ica_spatial.fig', output_folder) );
    saveas(gcf, sprintf('%s/ica_spatial.pdf', output_folder));
end
%%
figure('papersize', [7, 2]);
init_fig;
axes('position', [0.025, 0.05, 0.95, 4]);
imagesc([], [-0.3, 1.0]);
colorbar('southoutside');  set(gca, 'fontsize', 18);
axis off;
if export_fig
    saveas(gcf, sprintf('%s/colorbar_ica.fig', output_folder));
    saveas(gcf, sprintf('%s/colorbar_ica.pdf', output_folder));
    saveas(gcf, sprintf('%s/colorbar_ica.png', output_folder));
end
%%
figure('papersize', [7, 2]);
init_fig;
axes('position', [0.025, 0.05, 0.95, 4]);
imagesc([], [0, 1.0]);
colorbar('southoutside');  set(gca, 'fontsize', 18);
axis off;
if export_fig
    saveas(gcf, sprintf('%s/colorbar_cnmfe.fig', output_folder));
    saveas(gcf, sprintf('%s/colorbar_cnmfe.pdf', output_folder));
    saveas(gcf, sprintf('%s/colorbar_cnmfe.png', output_folder));
end

%% temporal traces
T = size(neuron.C, 2);
% temporal
figure('papersize', [9, 5]);
init_fig;
set(gcf, 'defaultAxesFontSize', 18);
col = colormap(cool);

C = neuron.C_raw(1:K_match,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_match:-1:1)');
axes('position', [0.08, 0.34, 0.89, 0.65]);
% plot((1:T)/neuron.Fs, C','linewidth', 1);
hold on;
ind_col = round(linspace(64, 1, size(C,1)));
for m=1:size(C,1)
    plot((1:T)/neuron.Fs, C(m,:), 'linewidth', 1, 'color', col(ind_col(m), :));
end
axis tight;
set(gca, 'ytick', 1:2:K_match);
set(gca, 'yticklabel', K_match:-2:1);
set(gca, 'xticklabel', []);
ylabel('Match');
box on ;
% xlim([0, 302]);

C = neuron.C_raw((K_match+1):K_cnmfe,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_cnmfe:-1:(K_match+1))');
axes('position', [0.08, 0.13, 0.89, 0.19]);
plot((1:T)/neuron.Fs, C', 'k', 'linewidth', 1);
axis tight;
set(gca, 'ytick', (K_match+1):1:K_cnmfe);
set(gca, 'yticklabel', K_cnmfe:-1:(K_match+1));
xlabel('Time (sec.)');
box on ;
ylabel('Other');
% xlim([0, 302]);

if export_fig
    saveas(gcf, sprintf('%s/cnmfe_temporal.fig', output_folder) );
    saveas(gcf, sprintf('%s/cnmfe_temporal.pdf', output_folder));
end

%% temporal traces of ICA neurons
% temporal
figure('papersize', [9, 5.0]);
init_fig;
set(gcf, 'defaultAxesFontSize', 18);

C = neuron_ica.C_raw(1:K_match,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_match:-1:1)');
axes('position', [0.08, 0.34, 0.89, 0.65]);
hold on;
for m=1:size(C,1)
    plot((1:T)/neuron.Fs, C(m,:), 'linewidth', 1, 'color', col(ind_col(m), :));
end
axis tight;
set(gca, 'ytick', 1:2:K_match);
set(gca, 'yticklabel', K_match:-2:1);
set(gca, 'xticklabel', []);
ylabel('Match');
box on ;

% xlim([0, 302]);
C = neuron_ica.C_raw((K_match+1):K_ica,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_ica:-1:(K_match+1))');
axes('position', [0.08, 0.13, 0.89, 0.19]);
plot((1:T)/neuron.Fs, C', 'k', 'linewidth', 1);
axis tight;
set(gca, 'ytick', (K_match+1):1:K_ica);
set(gca, 'yticklabel', K_ica:-1:(K_match+1));
xlabel('Time (sec.)');
box on ;
ylabel('Other');

% xlim([0, 302]);
if export_fig
    saveas(gcf, sprintf('%s/ica_temporal.fig', output_folder) );
    saveas(gcf, sprintf('%s/ica_temporal.pdf', output_folder));
end

%% temporal traces, cnmfe shorter.
% temporal
figure('papersize', [6, 4]);
init_fig;
set(gcf, 'defaultAxesFontSize', 18);
col = colormap(cool);

C = neuron.C_raw(1:K_match,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_match:-1:1)');
axes('position', [0.12, 0.27, 0.87, 0.722]);
% plot((1:T)/neuron.Fs, C','linewidth', 1);
hold on;
ind_col = round(linspace(64, 1, size(C,1)));
for m=1:size(C,1)
    plot((1:T)/neuron.Fs, C(m,:), 'linewidth', 1, 'color', col(ind_col(m), :));
end
axis tight;
set(gca, 'ytick', 1:2:K_match);
set(gca, 'yticklabel', K_match:-2:2);
set(gca, 'xticklabel', []);
ylabel('Match');
box on ;
xlim([0, 310]);

C = neuron.C_raw((K_match+1):K_cnmfe,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_cnmfe:-1:(K_match+1))');
axes('position', [0.12, 0.07, 0.87, 0.181]);
plot((1:T)/neuron.Fs, C', 'k', 'linewidth', 1);
axis tight;
set(gca, 'ytick', (K_match+1):1:K_cnmfe);
set(gca, 'yticklabel', K_cnmfe:-1:(K_match+1));
% xlabel('Time (sec.)');
set(gca, 'xticklabel', []);
box on ;
ylabel('Other');
xlim([0, 310]);

if export_fig
    saveas(gcf, sprintf('%s/cnmfe_temporal_v2.fig', output_folder) );
    saveas(gcf, sprintf('%s/cnmfe_temporal_v2.pdf', output_folder));
end

%% temporal traces of ICA neurons, shorter
% temporal
figure('papersize', [6,4]);
init_fig;
set(gcf, 'defaultAxesFontSize', 18);

C = neuron_ica.C_raw(1:K_match,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_match:-1:1)');
axes('position', [0.12, 0.27, 0.87, 0.722]);
hold on;
for m=1:size(C,1)
    plot((1:T)/neuron.Fs, C(m,:), 'linewidth', 1, 'color', col(ind_col(m), :));
end
axis tight;
set(gca, 'ytick', 1:2:K_match);
set(gca, 'yticklabel', K_match:-2:1);
set(gca, 'xticklabel', []);
% ylabel('Match');
box on ;

xlim([0, 310]);
C = neuron_ica.C_raw((K_match+1):K_ica,:);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (K_ica:-1:(K_match+1))');
axes('position', [0.12, 0.07, 0.87, 0.181]);
plot((1:T)/neuron.Fs, C', 'k', 'linewidth', 1);
axis tight;
set(gca, 'ytick', (K_match+1):2:K_ica);
set(gca, 'yticklabel', K_ica:-2:(K_match+1));
% xlabel('Time (sec.)');
box on ;
% ylabel('Other');
set(gca, 'xticklabel', []);
xlim([0, 310]);

axes('position', [0.12, 0.01, 0.87, 0.06]);
plot([280, 300], [1, 1], 'k', 'linewidth', 6);
axis off;
xlim([0, 310]);
if export_fig
    saveas(gcf, sprintf('%s/ica_temporal_v2.fig', output_folder) );
    saveas(gcf, sprintf('%s/ica_temporal_v2.pdf', output_folder));
end

%% save intervention results 
video_intervention; 


