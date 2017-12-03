%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra');
addpath(genpath('../../cbrewer')); 
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_BNST_subfigs'];

results_folder = sprintf('%sfig_BNST%s', results_folder, filesep); 
if ~exist(results_folder, 'dir')
    mkdir(results_folder); 
end  
results_file = [results_folder, 'fig_BNST_results.mat'];
cellsort_folder = [code_folder, 'CellSort'];
addpath(cellsort_folder);
demixed_video = [video_folder, filesep, 'BNST_decompose.avi'];
overlap_video = [video_folder, filesep, 'BNST_overlapping.avi'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file, 'Writable', true);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;


tshock = [690, 720, 780, 840, 930, 960, 1020, 1110, 1140, 1200, 1230];
pixel_size = 1000/726.67; %micron 
%% run CNMF-E or load results
nam = get_fullname('../../../data/CAMKII120_160317shock.mat');          % this demo data is very small, here we just use it as an example
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
nPCs = 200;
nICs = 150;
try
    neuron_ica = results_data.neuron_ica;
    A_ica_before_trim = results_data.A_ica_before_trim;
catch 
    run_ica;
    results_data.A_ica_before_trim = A_ica_before_trim;
    results_data.A_ica = A_ica;
    results_data.C_ica = C_ica;
    
    neuron_ica.compress_results();
    results_data.neuron_ica = neuron_ica;
    results_data.time_cost_ica = time_cost; 
end

%% compute SNRs 
K_cnmfe = size(neuron.C_raw, 1); 
K_ica = size(neuron_ica.C, 1); 
ind_frames = round(tshock(1)*neuron.Fs):round(tshock(end)*neuron.Fs); 
snr_ica = var(neuron_ica.C(:, ind_frames), 0, 2)./var(neuron_ica.C_raw(:, ind_frames)-neuron_ica.C(:, ind_frames), 0, 2); 
snr_cnmfe = var(neuron.C(:, ind_frames), 0, 2)./var(neuron.C_raw(:, ind_frames)-neuron.C(:, ind_frames), 0, 2); 

[~, srt] = sort(snr_cnmfe, 'descend'); 
neuron.orderROIs(srt); 
snr_cnmfe = snr_cnmfe(srt); 

[~, srt] = sort(snr_ica, 'descend'); 
neuron_ica.orderROIs(srt); 
A_ica_before_trim = A_ica_before_trim(:, srt); 
snr_ica = snr_ica(srt); 

%% find example neurons 
pixels = [35096, 10956, 12961, 9481, 19010, 10394, 17031, 35516, 15393, 25102, 13153, 38550]; 
ind_match = zeros(length(pixels), 1); 
A = neuron.A; 
A = bsxfun(@times, A, 1./max(A, [],1)); 
for m=1:length(pixels)
    temp = A(pixels(m),:); 
    [~, ind_match(m)] = max(temp); 
end
ind_match = sort(ind_match); 
%% 
% ind_match = [18, 14, 50, 27, 83, 73, 68, 10, 56, 100, 89, 74]; 
%% match neurons 
match_cnmfe = pair_neurons(neuron.A, neuron.C, neuron_ica.A, neuron_ica.C); 
ind_ica = match_cnmfe.ind_spatial(ind_match); 

%% compute correlation images 
d1 = neuron.options.d1; 
d2 = neuron.options.d2; 
create_correlation_images; 

%%  contours of the detected neurons, CNMF-E
neuron.Cn = Cn_filter;
neuron.PNR = PNR_filter;
if ~exist('Coor_cnmfe', 'var') || (length(Coor_cnmfe)~=size(neuron.A, 2))
    Coor_cnmfe = neuron.get_contours(0.8);
    neuron.Coor = Coor_cnmfe;
end
if ~exist('Coor_ica', 'var')|| (length(Coor_ica)~=size(neuron_ica.A, 2))
    Coor_ica = neuron_ica.get_contours(0.8);
    neuron_ica.Coor = Coor_ica;
end

%% show CNMF-E matched neurons
figure('papersize', [3,4]);
init_fig;
ctr = neuron.estCenter();
w = 2*neuron.options.gSiz + 1;
K_cnmfe = size(neuron.A,2);
K_ica = size(neuron_ica.A, 2);
K_match = length(ind_match);
dw = 1/3;
dh = 1/4;
for m=1:12
  ii = mod(m, 3);
    jj = ceil(m/3);
    if ii==0
        ii = 3;  
    end 
    axes('position', [(ii-1)*dw+0.002, 1-jj*dh+0.001, dw*0.95, dh*0.95]);
    
    
    r0 = round(ctr(ind_match(m),1));
    c0 = round(ctr(ind_match(m),2));
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
    
    img = neuron.reshape(neuron.A(:, ind_match(m)), 2);
    temp = img(r0+(-w:w), c0+(-w:w)); 
    imagesc(temp/max(temp(:)), [0,1]);
    axis equal off; hold on;
    text(1, 8, num2str(m), 'color', 'w','fontweight', 'bold','fontsize', 22);
end
if export_fig
    saveas(gcf, sprintf('%s/cnmfe_match_spatial.fig', output_folder));
    saveas(gcf, sprintf('%s/cnmfe_match_spatial.pdf', output_folder));
end
% temporal traces
figure('papersize', [6.8, 1.8]);
init_fig;
set(gcf, 'defaultAxesFontSize', 11);
col = colormap(cool); 
ind = round(650*neuron.Fs):round(1250*neuron.Fs);
C = neuron.C_raw(ind_match,ind);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (12:-1:1)');
axes('position', [0.05, 0.19, 0.9, 0.8]); hold on; 
for m=1:12
plot((0:(size(C,2)-1))/neuron.Fs, C(m,:), 'color', col(64-4*m, :),'linewidth', 1);
end 
axis tight;
set(gca, 'ytick', 1:2:12);
set(gca, 'yticklabel', 12:-2:1);
% set(gca, 'xticklabel', []);
% xlabel('Time (sec.)');
box on; 
set(gca, 'xtick', []); 
xlim([0, 600]);
ylim([0.5, 13.5]);
for m=1:length(tshock)
    plot(ones(2,1)*tshock(m)-650, get(gca, 'ylim'), '-.k', 'linewidth', 0.5);
end
if export_fig
    saveas(gcf, sprintf('%s/cnmfe_match_temporal.fig', output_folder));
    saveas(gcf, sprintf('%s/cnmfe_match_temporal.pdf', output_folder));
end

%% show ICA matched neurons
figure('papersize', [3,4]);
init_fig;
ctr = neuron.estCenter();
w = 2*neuron.options.gSiz + 1;
K_cnmfe = size(neuron.A,2);
K_ica = size(neuron_ica.A, 2);
K_match = length(ind_match);
dw = 1/3;
dh = 1/4;
for m=1:12
    ii = mod(m, 3);
    jj = ceil(m/3);
    if ii==0
        ii = 3;  
    end 
    axes('position', [(ii-1)*dw+0.002, 1-jj*dh+0.001, dw*0.95, dh*0.95]);
    
    
    r0 = round(ctr(ind_match(m),1));
    c0 = round(ctr(ind_match(m),2));
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
    
    img = neuron.reshape(A_ica_before_trim(:, ind_ica(m)), 2);
    temp = img(r0+(-w:w), c0+(-w:w)); 
    imagesc(temp/max(temp(:)), [-0.3,1]);
    axis equal off; hold on;
    text(1, 8, num2str(m), 'color', 'w','fontweight', 'bold', 'fontsize', 22);
end
if export_fig
    saveas(gcf, sprintf('%s/ica_match_spatial.fig', output_folder));
    saveas(gcf, sprintf('%s/ica_match_spatial.pdf', output_folder));
end
% temporal traces
figure('papersize', [6.8, 1.8]);
init_fig;
set(gcf, 'defaultAxesFontSize', 11);

ind = round(650*neuron.Fs):round(1250*neuron.Fs);
C = neuron_ica.C_raw(ind_ica,ind);
C = bsxfun(@times, C, 1.3./max(C,[],2));
C = bsxfun(@plus, C, (12:-1:1)');
axes('position', [0.05, 0.19, 0.9, 0.8]); hold on; 
for m=1:12
plot((0:(size(C,2)-1))/neuron.Fs, C(m,:), 'color', col(64-4*m, :),'linewidth', 1);
end 
axis tight;
set(gca, 'ytick', 1:2:12);
set(gca, 'yticklabel', 12:-2:1);
% set(gca, 'xticklabel', []);
% xlabel('Time (sec.)');
box on ;
set(gca, 'xtick', []); 
xlim([0, 600]);
ylim([0.5, 13.5]);

for m=1:length(tshock)
    plot(ones(2,1)*tshock(m)-650, get(gca, 'ylim'), '-.k', 'linewidth', 0.5);
end

axes('position', [0.05, 0.08, 0.9, 0.1]); 
plot([591, 600], [1, 1], 'k', 'linewidth', 3); 
xlim([0, 600]); 
axis off; 
if export_fig
    saveas(gcf, sprintf('%s/ica_match_temporal.fig', output_folder));
    saveas(gcf, sprintf('%s/ica_match_temporal.pdf', output_folder));
end

%% 
figure('papersize', [0.9, 4]); 
init_fig; 
axes('position', [-2.08, 0.05, 2.8, 0.9]); 
set(gca, 'fontsize', 15); 
imagesc([], [-0,1]); 
colorbar;
axis off; 
if export_fig
    saveas(gcf, sprintf('%s/cnmfe_colorbar.fig', output_folder));
    saveas(gcf, sprintf('%s/cnmfe_colorbar.png', output_folder));
    saveas(gcf, sprintf('%s/cnmfe_colorbar.pdf', output_folder));
end

figure('papersize', [0.9, 4]); 
init_fig; 
axes('position', [-2.08, 0.05, 2.8, 0.9]); 
set(gca, 'fontsize', 15); 
imagesc([], [-0.3,1]); 
colorbar;
axis off; 
if export_fig
    saveas(gcf, sprintf('%s/ica_colorbar.fig', output_folder));
    saveas(gcf, sprintf('%s/ica_colorbar.png', output_folder));
    saveas(gcf, sprintf('%s/ica_colorbar.pdf', output_folder));
end
%% plot contours
figure('papersize', [d2+50, d1+20]/max(d1,d2)*5);
init_fig;

imagesc(Cn_filter, [0.5, 1]); colormap gray; hold on; colorbar;
axis equal off tight;
set(gca,'position',[0.01 0.08 .8 0.9],'units','normalized'); 
for m=1:length(ind_match)
    temp = Coor_cnmfe{ind_match(m)};
    fill(temp(1,:), temp(2,:), 'r', 'facealpha', 0.5,'edgecolor', 'w');
    [~, ind] = max(temp(1,:)); %randi(length(temp), 1);
    text(temp(1,ind), temp(2,ind), num2str(m), 'color', 'r', ...
        'fontsize', 14,'fontweight', 'bold');
end

for m=1:length(ind_match)
    temp = Coor_ica{ind_ica(m)};
    fill(temp(1,:), temp(2,:), 'g', 'facealpha', 0.5,'edgecolor', 'w');
    [~, ind] = min(temp(1,:)); %randi(length(temp), 1);
    text(temp(1,ind)-7, temp(2,ind), num2str(m), 'color', 'g', ...
        'fontsize', 14, 'fontweight', 'bold');
end
a = plot([220, 260]/pixel_size, 190*[1,1], 'c', 'linewidth', 10); 
b = text(220/pixel_size, 180, '40 um', 'color', 'c', 'fontsize', 18, 'fontweight', 'bold');

axes('position',[0.15, 0.02, .8, 0.07],'units','normalized'); hold on; 
temp1 = fill(ones(3,1), ones(3,1), 'r', 'facealpha', 0.5, 'edgecolor', 'w'); 
temp2 = fill(ones(3,1), ones(3,1), 'g', 'facealpha', 0.5, 'edgecolor', 'w'); 
legend('CNMF-E', 'PCA/ICA', 'orientation', 'horizental'); 
axis equal off tight; box off; legend boxoff; 
if export_fig
    saveas(gcf, sprintf('%s/match_contours.fig', output_folder));
    saveas(gcf, sprintf('%s/match_contours.pdf', output_folder));
end


%% show median activity after shocks 
ind = bsxfun(@plus, (-10:40)', round(neuron.Fs*tshock))-1;
xx = (-10:40)/neuron.Fs;
yy = (1:12);
cmap_red = cbrewer('seq', 'Reds', 64); 
cmap_green = cbrewer('seq', 'Greens', 64); 
h1 = 0.2; 
v_alpha = 0.4; 
peak_cnmfe = zeros(12, 1); 
mad_cnmfe = zeros(12,1); 
peak_ica = zeros(12,1); 
mad_ica = zeros(12,1); 

post_cnmfe = zeros(12, 11); 
pre_cnmfe = zeros(12,11); 
post_ica = zeros(12,11); 
pre_ica = zeros(12,11); 

for m=1:length(ind_match)
    if m==1
        figure('papersize', [2, 7.2]);
    else
        figure('papersize', [1.2, 7.2]);
    end
    init_fig;
    set(gcf, 'defaultAxesFontSize', 20);
    cc = neuron.C_raw(ind_match(m), :);
    cc_ica = neuron_ica.C_raw(ind_ica(m), :);
    temp = cc(ind);
    temp_ica = cc_ica(ind);
    %  CNMF-E 
    imagesc(xx, yy, temp', quantile(temp(:), [0.01, 1]));
    hold on; axis tight xy ;
    set(gca, 'xtick', 0:5:15);
    set(gca, 'xticklabel', []);
    hold on;
    plot([0,0], get(gca, 'ylim'), '-.k');
    tmp_xlim = get(gca, 'xlim');
    set(gca, 'xtick', 0:5:15);
    set(gca, 'xticklabel', []);
        colormap(gca, cmap_red); 
    if m==1
        ylabel('Trial #');
        set(gca, 'position', [0.402, 0.93-h1, 0.588, h1]);
        set(gca, 'ytick', 3.5:3:9.5);
        set(gca, 'yticklabel', 3:3:9);
        title(sprintf('Cell # %d', m));
        axes('position', [0.402, 0.5, 0.588, 0.43-h1]);
    else
        set(gca, 'position', [0.01, 0.93-h1, 0.98, h1]);
        title(sprintf('# %d', m));
        axes('position', [0.01, 0.5, 0.98, 0.43-h1]);
    end
    
    ymean = median(temp,2);
    ste = median(abs(bsxfun(@minus, temp, ymean)), 2); %norminv(0.975)*std(temp, 0, 2)/sqrt(size(temp,2));
    fill([xx'; flipud(xx')], [ymean-ste; flipud(ymean+ste)], 'r', 'facealpha', ...
        v_alpha, 'edgecolor', 'none'); hold on;
    plot(xx, ymean, 'color', cmap_red(end,:), 'linewidth', 2);
    axis  off;     ylim(mean(ymean)+(max(ymean+ste)-mean(ymean))*[-0.8,1.02]); 
    plot([0,0], get(gca, 'ylim'), '-.k');
    axis off;
    xlim(tmp_xlim);
     post_cnmfe(m, :) = mean(temp(12:16, :)); 
    pre_cnmfe(m, :) = mean(temp(6:10, :)); 
    
    peak_cnmfe(m) = max(ymean)-mean(ymean(1:10)); 
    mad_cnmfe(m) = mean(ste); 
    
    % PCA/ICA 
    axes;    imagesc(xx, yy, temp_ica', quantile(temp_ica(:), [0.01, 1]));
      hold on; axis tight xy ;
    set(gca, 'xtick', 0:5:15);
    set(gca, 'xticklabel', []);
    hold on;
    plot([0,0], get(gca, 'ylim'), '-.k');
    tmp_xlim = get(gca, 'xlim');
    set(gca, 'xtick', 0:5:15);
    set(gca, 'xticklabel', []);
    colormap(gca, cmap_green); 
    if m==1
        set(gca, 'position', [0.402, 0.51-h1, 0.588, h1]);
                    set(gca, 'ytick', 3.5:3:9.5);     
                    set(gca, 'yticklabel', 3:3:9); 
        ylabel('Trial #');
        set(gca, 'xtick', 0:5:15);
        set(gca, 'xticklabel', []);
        axes('position', [0.402, 0.08, 0.588, 0.43-h1]);
    else
        set(gca, 'position', [0.01, 0.51-h1, 0.98, h1]);
        set(gca, 'ytick', 3.5:3:9.5);
        set(gca, 'yticklabel', []); 
        set(gca, 'xtick', 0:5:15);
        set(gca, 'xticklabel', []);
        axes('position', [0.01, 0.08, 0.98, 0.43-h1]);
    end
    %     xlabel('Time (sec.)');
    ymean = median(temp_ica,2);
    ste = median(abs(bsxfun(@minus, temp_ica, ymean)), 2); %norminv(0.975)*std(temp, 0, 2)/sqrt(size(temp,2));
    fill([xx'; flipud(xx')], [ymean-ste; flipud(ymean+ste)], 'g', 'facealpha', ...
        v_alpha, 'edgecolor', 'none');  hold on;
    plot(xx, ymean, 'color', cmap_green(end, :), 'linewidth', 2);
    axis  off;
    ylim(mean(ymean)+(max(ymean+ste)-mean(ymean))*[-0.8, 1.02]); 
    plot([0,0], get(gca, 'ylim'), '-.k');
    xlim(tmp_xlim);
    
     post_ica(m, :) = mean(temp_ica(12:16, :)); 
    pre_ica(m, :) = mean(temp_ica(6:10, :)); 
     peak_ica(m) = max(ymean)-mean(ymean(1:10)); 
    mad_ica(m) = mean(ste); 
    if m==12
        axes('position', [0.01, 0.03, 0.98, 0.03]);
        plot(xx(end)+[-10,-1], [0, 0], 'k', 'linewidth', 8);
        xlim(tmp_xlim); axis off; 
    end
    
    if export_fig
        saveas(gcf, sprintf('%s/shock_neuron_%d.fig', output_folder, m));
        saveas(gcf, sprintf('%s/shock_neuron_%d.pdf', output_folder, m));
    end
end

%%
figure('papersize', [4.8, 7]);
init_fig; 
set(gca, 'fontsize', 25); 
hold on; 
col = colormap(cool); 
for m=1:12
plot(peak_ica(m)/mad_ica(m), peak_cnmfe(m)/mad_cnmfe(m), 'ow', ...
    'markerfacecolor', col(64-4*m, :), 'markersize', 15 ); 
text(peak_ica(m)/mad_ica(m)+0.1, peak_cnmfe(m)/mad_cnmfe(m)+0.2, num2str(m), 'fontsize', 25);
end
box on; axis tight; 
plot([0, 5.5], [0, 5.5], 'r'); 
ylim([0, 11]);
set(gca, 'xtick', 0:2:4); 
xlabel('PCA/ICA'); 
ylabel('CNMF-E'); 
if export_fig
    saveas(gcf, sprintf('%s/shock_peak.fig', output_folder));
    saveas(gcf, sprintf('%s/shock_peak.pdf', output_folder));
end
%%
% diff_cnmfe = post_cnmfe-pre_cnmfe;
% diff_ica = post_ica-pre_ica; 
% med_cnmfe = median(diff_cnmfe,2); 
% med_ica = median(diff_ica,2); 
% mad_cnmfe = median(abs(bsxfun(@minus, diff_cnmfe, med_cnmfe)), 2); 
% mad_ica = median(abs(bsxfun(@minus, diff_ica, med_ica)), 2); 
% figure('papersize', [6,6]);
% init_fig; 
% hold on; 
% col = colormap(cool); 
% for m=1:12
% plot(med_ica(m)/mad_ica(m), med_cnmfe(m)/mad_cnmfe(m), 'ow', ...
%     'markerfacecolor', col(64-4*m, :), 'markersize', 10 ); 
% text(med_ica(m)/mad_ica(m)-0.15, med_cnmfe(m)/mad_cnmfe(m)+0.1, num2str(m), 'fontsize', 16);
% end
% box on; axis tight; 
% plot([-1, 2.5], [-1, 2.5], 'r'); 
% if export_fig
%     saveas(gcf, sprintf('%s/shock_diff.fig', output_folder));
%     saveas(gcf, sprintf('%s/shock_diff.pdf', output_folder));
% end

%% save video
kt = 1;     % play one frame in every kt frames
save_avi = true;
t_begin = round(650*neuron.Fs); 
t_end = round(1250*neuron.Fs); 
ind_shock = round(tshock*neuron.Fs); 
k_shock = 0; 
next_shock = ind_shock(1); 
last_shock = -inf; 
center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2));
range_res = [-1,1]*center_ac;
range_ac = center_ac+range_res;
multi_factor = 3;
center_Y = min(Y(:)) + multi_factor*center_ac;
range_Y = center_Y + range_res*multi_factor;

avi_filename = demixed_video; %sprintf('../../Videos/footshock_decompose.avi'); 
cnmfe_save_video; 


% 
% %% demixing video 
% file_nm = sprintf('../../Videos/footshock_demix.avi'); 
% neuron.playAC(file_nm, [], [t_begin, t_end]); 
% 
% 



