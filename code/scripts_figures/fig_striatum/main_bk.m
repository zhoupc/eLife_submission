%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_Vessel_subfigs'];
results_file = [results_folder, 'fig_striatum.mat'];
cellsort_folder = [code_folder, 'CellSort'];
addpath(cellsort_folder);
demixed_video = [video_folder, filesep, 'striatum_decompose.avi'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file, 'Writable', true);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% run CNMF-E or load results
nam = get_fullname('../../../data/blood_vessel_10Hz.mat');          % this demo data is very small, here we just use it as an example
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
if ~exist('Y', 'var')
    Y = neuron.load_patch_data();
end
Y = neuron.reshape(Y, 1);
Ybg = neuron.reconstruct_background();
Ybg = neuron.reshape(Ybg, 1);

%% run PCA/ICA analysis or load results
nPCs = 3000;
nICs = 700;
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
end

%% compute SNRs and sort CNMF-E and PCA/ICA neurons according to its SNR

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

%% match neurons and reorder the sequence
% match_cnmfe = pair_neurons(neuron.A, neuron.C, neuron_ica.A, neuron_ica.C_raw);
% ids_match = find(~isnan(match_cnmfe.ind_max));
% ids_ica = match_cnmfe.ind_max(ids_match);
% K_match = length(ids_match);
%
% ind = 1:K_cnmfe;
% ind(ids_match) = [];
% srt = [ids_match, ind];
% neuron.orderROIs(srt);
% snr_cnmfe = snr_cnmfe(srt);
%
% ind = 1:K_ica;
% ind(ids_ica) = [];
% srt = [ids_ica, ind];
% neuron_ica.orderROIs(srt);
% snr_ica = snr_ica(srt);
% A_ica_before_trim = A_ica_before_trim(:, srt);

%%
match_ica = pair_neurons(neuron_ica.A, neuron_ica.C_raw, neuron.A, neuron.C);

ids_ica = find(~isnan(match_ica.ind_max));
ids_match = match_ica.ind_max(ids_ica);
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

%% show one frame and one neuron
Ymin = 1500;
Ymax = 3000;
Yacmin = 10;
Yacmax = 200;
d1 = neuron.options.d1;
d2 = neuron.options.d2;
[~, example_neuron] = max(neuron.A(82*256+172, :));
ind_frame = 1082;
% [~, ind_frame] = max(neuron.C(example_neuron, :));
figure;
% Cn_res = correlation_image(Ysignal-neuron.A*neuron.C, 8, d1, d2);
coor = plot_contours(neuron.A(:,example_neuron), randn(d1,d2), 0.9, 0);
close;
% create a mask  for selecting ROI
x = round(coor{1}(1,5:end));
y = round(coor{1}(2,5:end));
mask = false(d1,d2);
for ii =1:length(x)
    for jj=1:length(x)
        tmp_x = 0: (x(jj)-x(ii));
        if length(tmp_x)<2
            continue;
        end
        tmp_y = y(ii)+(y(jj)-y(ii))/(x(jj)-x(ii))*tmp_x;
        mask(sub2ind([d1,d2], round(tmp_y), tmp_x+x(ii))) = true;
    end
end

y_raw = mean(Y(mask(:), 2001:3000));
y_bg = mean(Ybg(mask(:), 2001:3000));
baseline = mean2(Ybg(mask(:), :));
t = (1:1000)/neuron.Fs;

figure('papersize', [6, 4]);
init_fig;
set(gcf,  'defaultAxesFontSize',14);

dy = 350;
axes('position', [0.01, 0.01, 0.98, 0.98]);

plot(t, y_raw-baseline+2*dy, 'k', 'linewidth', 1);
hold on;
plot(t, y_bg-baseline+dy, 'b', 'linewidth', 1);
plot(t, y_raw-y_bg, 'r', 'linewidth', 1);
plot([t(end)-8, t(end)+2], [-40, -40], 'k', 'linewidth', 6);
% text(t(end)-9,-72, '10 sec.', 'fontsize', 12);
% temp = text(t(end)+7, -55, '100 (a.u.)', 'fontsize', 12');
% set(temp, 'rotation', 90);
axis off;
% text(60, 2*dy+180, sprintf('variance: %.2f', var(y_raw)), 'fontsize', 12)
% text(60, dy+180, sprintf('variance: %.2f', var(y_bg)),'color', 'b', 'fontsize', 12)
% text(60, 180, sprintf('variance: %.2f', var(y_raw-y_bg)),'color', 'r', 'fontsize', 12)
text(1, 2*dy+180, 'raw data', 'fontsize', 20);
text(1, dy+180, 'estimated background', 'color', 'b', 'fontsize', 20);
text(1, 180, 'background-subtracted data','color', 'r', 'fontsize', 20);
xmin = -3;
xmax = t(end)+5;
ymin = -100;
ymax = 2*dy+220;
plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k');
ylim([ymin, ymax]+[-1, 1]);
xlim([xmin, xmax]+[-1,1]);
%
% set(gca, 'xticklabel', []);
% ylim([-110,290]);
% legend('estimated background')
% ylabel('Fluo.');
% text(5, 230, sprintf('variance: %.2f', var(y_bg)), 'fontsize', 12)
%
% axes('position', [0.13, 0.15, 0.82, 0.24]);
% plot(t, y_raw-y_bg, 'r', 'linewidth', 1);
% xlabel('Time (sec.)');
% ylim([-60,340]); hold on;
% legend('background-subtracted data');
% ylabel('Fluo.');
% text(5, 280, sprintf('variance: %.2f', var(y_raw-y_bg)), 'fontsize', 12)

% axes('position', [0.13, 0.75, 0.82, 0.24]);
%
% plot(t, y_raw-baseline, 'k');
% set(gca, 'xticklabel', [], 'linewidth', 1);
% ylim([-110,290]);
% legend('raw data');
% ylabel('Fluo.');
% text(5, 230, sprintf('variance: %.2f', var(y_raw)), 'fontsize', 12)
%
% axes('position', [0.13, 0.45, 0.82, 0.24]);
% plot(t, y_bg-baseline, 'b', 'linewidth', 1);
% set(gca, 'xticklabel', []);
% ylim([-110,290]);
% legend('estimated background')
% ylabel('Fluo.');
% text(5, 230, sprintf('variance: %.2f', var(y_bg)), 'fontsize', 12)
%
% axes('position', [0.13, 0.15, 0.82, 0.24]);
% plot(t, y_raw-y_bg, 'r', 'linewidth', 1);
% xlabel('Time (sec.)');
% ylim([-60,340]); hold on;
% legend('background-subtracted data');
% ylabel('Fluo.');
% text(5, 280, sprintf('variance: %.2f', var(y_raw-y_bg)), 'fontsize', 12)

if export_fig
    saveas(gcf, sprintf('%s/example_roi_traces.fig', output_folder) );
    saveas(gcf, sprintf('%s/example_roi_traces.pdf', output_folder));
end

%% show decomposed images for one example frame
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
Y = neuron.reshape(Y, 1);
Ybg = neuron.reshape(Ybg, 1);

set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
neuron.image(Y(:, ind_frame), [Ymin, Ymax]);
colorbar;
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
plot([1,1,128, 128, 1], [129, 256, 256, 129, 129], '-.m');
axis equal off tight;
if export_fig
    saveas(gcf, sprintf('%s/example_frame.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end

figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ybg(:, ind_frame), [Ymin, Ymax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized');
colorbar;
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_bg.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_bg.pdf', output_folder));
end

b0 = mean(Ybg,2);
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ybg(:, ind_frame)-b0); %, [-50, 150]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized');
colorbar;
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_bg_fluc.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_bg_fluc.pdf', output_folder));
end


figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(b0, [Ymin, Ymax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized');
colorbar;
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_bg_constant.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_bg_constant.pdf', output_folder));
end


figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(neuron.A*neuron.C(:, ind_frame), [Yacmin, Yacmax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized');
colorbar;
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_ac.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_ac.pdf', output_folder));
end

figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(double(Y(:, ind_frame))-Ybg(:, ind_frame)- neuron.A*neuron.C(:, ind_frame), [-0.5, 0.5]*Yacmax);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2);
axis equal off tight; colorbar;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_res.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_res.pdf', output_folder));
end

%% variance reduced
try
    var_0 = results_data.var0;
    var_nobg = results_data.var_nobg;
    var_final = results_data.var_nobg;
    var_no_signal = results_data.var_no_signal;
catch
    var_0 = var(double(Y), [], 2);  % variance in the data
    var_nobg = var(double(Y)-Ybg,[],2);
    var_final = var(double(Y)-Ybg-neuron.A*neuron.C, [], 2);
    var_no_signal = var(double(Y)-neuron.A*neuron.C, [],2);
    
    results_data.var0 = var_0;
    results_data.var_nobg = var_nobg;
    results_data.var_nobg = var_final;
    results_data.var_no_signal = var_no_signal;
end
% [u,s,v] = svdsecon(bsxfun(@minus, Ybg, mean(Ybg, 2)), 20);
% var_rank1 = var(double(Y)-u(:,1)*s(1)*v(:,1)', [],2);
% var_rank2 = var(double(Y)-u(:,1:2)*s(1:2,1:2)*v(:,1:2)', [], 2);
% var_rank3 = var(double(Y)-u(:,1:3)*s(1:3,1:3)*v(:,1:3)', [], 2);
% var_rank5 = var(double(Y)-u(:,1:5)*s(1:5,1:5)*v(:,1:5)', [], 2);
% var_rank10 = var(double(Y)-u(:,1:10)*s(1:10,1:10)*v(:,1:10)', [], 2);

%%
figure('papersize', [6, 4]);
init_fig;
set(gcf,  'defaultAxesFontSize',16);
set(gca, 'position', [0.11, 0.16, 0.87, 0.76]);
v_alpha = 0.5;
hold on ;
dbin = 0.005;
bins = 0:dbin:1;
% variance explained by background
[count, bins] = hist(1-var_nobg./var_0, bins);
tmp_p = count/sum(count)/dbin;
fill([bins, bins(end), bins(1)], [tmp_p, dbin, dbin],'b', 'facealpha', v_alpha, 'edgecolor', 'none');

% variance explained by the 1st PC of the background
% [count, bins] = hist(1-var_rank1./var_0, 50);
% fill([bins, bins(end), bins(1)], [count, 1,1],'m', 'facealpha', v_alpha, 'edgecolor', 'none');

% variance explained by 2nd-10th PCs of the background
% [count, bins] = hist((var_rank1-var_rank3)./var_0, 50);
% fill([bins, bins(end), bins(1)], [count, 1,1],'m', 'facealpha', v_alpha, 'edgecolor', 'none');

% variance explained by all other background components except the first 1
% [count, bins] = hist((var_rank1-var_nobg)./var_0, 50);
% fill([bins, bins(end), bins(1)], [count, 1,1],'g', 'facealpha', v_alpha, 'edgecolor', 'none');

% variance explained by neural signal
[count, bins] = hist((var_nobg-var_final)./var_0, bins);
tmp_p = count/sum(count)/dbin;
fill([bins, bins(end), bins(1)], [tmp_p, dbin, dbin], 'r', 'facealpha', v_alpha, 'edgecolor', 'none');

% variance in the residual
[count, bins] = hist(var_final./var_0, bins);
tmp_p = count/sum(count)/dbin;
fill([bins, bins(end), bins(1)], [tmp_p, dbin, dbin],'g', 'facealpha', v_alpha, 'edgecolor', 'none');

box on;
xlim([0, 1]);
xlabel('Relative variance')
ylabel('Probability density');
temp = breakxaxis([0.2, 0.6]);
% legend('BG', '1st PC', 'BG without the 1st PC', 'denoised neural signal', 'residual');
axes(temp.rightAxes);
legend('background', 'denoised neural signal', 'residual');
axis off;
if export_fig
    saveas(gcf, sprintf('%s/variance_explained.fig', output_folder) );
    saveas(gcf, sprintf('%s/variance_explained.pdf', output_folder));
end

%% compute correlation images
create_correlation_images;

%%  contours of the detected neurons, CNMF-E
neuron.Cn = Cn_filter;
neuron.PNR = PNR_filter;
if ~exist('Coor_cnmfe', 'var')
    Coor_cnmfe = neuron.get_contours(0.6);
    neuron.Coor = Coor_cnmfe;
end
if ~exist('Coor_ica', 'var')
    Coor_ica = neuron_ica.get_contours(0.8);
    neuron_ica.Coor = Coor_ica;
end


%% missed neurons by PCA/ICA
ind_ica_miss = 1:K_cnmfe;
ind_ica_miss(match_ica.ind_spatial) = [];

%% contour plots of CNMF-E and PCA/ICA results.
ctr = neuron.estCenter();
figure('papersize', [6, d2/d1*5]);
init_fig;
[~, srt] = sort(snr_cnmfe, 'descend');
plot_contours(neuron.A(:, srt), Cn_filter, 0.6, 0, [], Coor_cnmfe(srt), 3);
set(gca, 'position', [0.01, 0.02, 0.9, 0.96]);

for m=ind_ica_miss
    %     temp = Coor_cnmfe{ind_miss(m)};
    if ~(ctr(m,1)>129 && ctr(m,2)<128 )
        continue;
    end
    
    temp = Coor_cnmfe{m};
    temp = medfilt1(temp')';
    h = fill(temp(1,3:end), temp(2,3:end), 'g', 'facealpha',0.3, 'edgecolor', 'none');
    uistack(h, 'bottom');
end

temp = get(gca, 'children');
for m=1:length(temp)
    h = temp(m);
    if strcmpi(h.Type, 'Image')
        uistack(h, 'bottom');
    end
end
plot([1,1,128, 128, 1], [129, 256, 256, 129, 129], 'm');
colormap gray; axis  off equal;
axis([0, 128, 129, 256]); colorbar;
if export_fig
    saveas(gcf, sprintf('%s/contours_cnmfe.fig', output_folder) );
    saveas(gcf, sprintf('%s/contours_cnmfe.pdf', output_folder));
    saveas(gcf, sprintf('%s/contours_cnmfe.png', output_folder));
    saveas(gcf, sprintf('%s/contours_cnmfe.fig', output_folder));
end

%%
figure('papersize', [6, d2/d1*5]);
ctr = neuron.estCenter();
init_fig;
[~, srt] = sort(snr_ica, 'descend');
plot_contours(neuron_ica.A(:, srt), Cn_filter, 0.6, 0, [], Coor_ica(srt), 3);
set(gca, 'position', [0.01, 0.02, 0.9, 0.96]);

for m=ind_ica_miss
    %     temp = Coor_cnmfe{ind_miss(m)};
    if ~(ctr(m,1)>129 && ctr(m,2)<128 )
        continue;
    end
    temp = Coor_cnmfe{m};
    temp = medfilt1(temp')';
    h = fill(temp(1,3:end), temp(2,3:end), 'g', 'facealpha',0.3, 'edgecolor', 'none');
    uistack(h, 'bottom');
end

temp = get(gca, 'children');
for m=1:length(temp)
    h = temp(m);
    if strcmpi(h.Type, 'Image')
        uistack(h, 'bottom');
    end
end
% for m=1:length(Coor_ica)
%     temp = Coor_ica{m};
%     plot(temp(1, :), temp(2,:), '-', 'color', [0,1, 0],  'linewidth', 2);
% end
colormap gray; axis tight off equal;
plot([1,1,128, 128, 1], [129, 256, 256, 129, 129], 'm');
axis([0, 128, 129, 256]);
colorbar;
if export_fig
    saveas(gcf, sprintf('%s/contours_ica.fig', output_folder) );
    saveas(gcf, sprintf('%s/contours_ica.pdf', output_folder));
    saveas(gcf, sprintf('%s/contours_ica.png', output_folder));
    saveas(gcf, sprintf('%s/contours_ica.fig', output_folder));
end

%% example of missed neurons, PCA/ICA
T = size(neuron.C, 2);
figure('papersize', [13,1]);
init_fig;
ctr = neuron.estCenter();

ind_miss = zeros(14, 1);
figure;
k = 0;
for m=(K_match+1):K_cnmfe
    subplot(221);
    neuron.image(neuron.A(:,m));
    axis equal off tight;
    subplot(2, 2, 3:4);
    plot(neuron.C_raw(m, :));
    xlim([0, 3000]);
    x = input('use this neuron as an example (y/n): ', 's');
    if strcmpi(x, 'y')
        k = k+1;
        ind_miss(k) = m;
        if k==14
            break;
        end
    end
end
% ind_miss = [389, 391, 394, 395, 397, 406, 408, 430,...
%     464, 473,491,476,496,536];

figure('papersize', [13,1]);
init_fig;
K_show = 14;
for m=1:K_show
    %     if m<5
    %        axes('position', [(m-1)*0.2+0.001, 0.501, 0.198, 0.498]);
    %     else
    %        axes('position', [(m-6)*0.2+0.001, 0.001, 0.198, 0.498]);
    %     end
    axes('position', [(m-1)*1/K_show+0.001, 0.001, 1/K_show-0.001, 0.98]);
    r0 = round(ctr(ind_miss(m), 1));
    c0 = round(ctr(ind_miss(m), 2));
    img = neuron.reshape(neuron.A(:, ind_miss(m)),2);
    if r0<23
        r0 = 23;
    elseif r0>=d1-18
        r0 = d1-18;
    end
    if c0<21
        c0 = 21;
    elseif c0>=d2-20
        c0 = d2-20;
    end
    imagesc(img(r0+(-22:18), c0+(-20:20)));
    
    axis equal off tight;
    text(1,6, num2str(m), 'color', 'w', 'fontweight', 'bold', 'fontsize', 15);
end
if export_fig
    saveas(gcf, sprintf('%s/ica_missed_spatial.fig', output_folder) );
    saveas(gcf, sprintf('%s/ica_missed_spatial.pdf', output_folder));
end

% temporal
figure('papersize', [14, 4.2]);
init_fig;
set(gcf, 'defaultAxesFontSize', 16);
col = colormap(cool);
C = neuron.C_raw(ind_miss(1:K_show),:);
axes('position', [0.0, 0.0, 0.95, 1]); hold on;
for m=1:K_show
    y = C(m,:);
    plot((1:T)/neuron.Fs, 1.2*y/max(y)+K_show-m,'linewidth', 1, 'color', [1,1,1]*0.3);
    text(-10-2*(m>=10), (K_show-m), num2str(m), 'fontsize', 20);
end
plot([290,300], [-0.5, -0.5], 'k', 'linewidth', 5);
% text(289, -1, '10 sec', 'fontsize', 14);
axis off;
axis([-20,300,-1.5, K_show]);
if export_fig
    saveas(gcf, sprintf('%s/ica_missed_temporal.fig', output_folder) );
    saveas(gcf, sprintf('%s/ica_missed_temporal.pdf', output_folder));
end

%% example of matched neurons
% spatial colorbar
% CNMF-E
figure('papersize', [2, 0.4]);
init_fig;
axes('position', [0.025, 0.1, 0.95, 0.85]);
colorbar('southoutside');  set(gca, 'fontsize', 20);
axis off;
saveas(gcf, sprintf('%s/colorbar_cnmfe.fig', output_folder));
saveas(gcf, sprintf('%s/colorbar_cnmfe.pdf', output_folder));

% PCA/ICA
figure('papersize', [2, 0.3]);
init_fig;
axes('position', [0.05, 0.05, 0.9, 4]);
imagesc([], [-0.2, 1.0]);
colorbar('southoutside');  set(gca, 'fontsize', 18);
axis off;
saveas(gcf, sprintf('%s/colorbar_ica.fig', output_folder));
saveas(gcf, sprintf('%s/colorbar_ica.pdf', output_folder));
% %%
% for m=351:K_match
%     subplot(321);
%     neuron.image(neuron.A(:, m));
%     axis equal off tight;
%     xlim(ctr(m, 2)+[-20, 20]);
%     ylim(ctr(m,1)+[-20, 20]);
%     subplot(322);
%     neuron.image(A_ica_before_trim(:, ));
%     axis equal off tight;
%
%     xlim(ctr(m, 2)+[-20, 20]);
%     ylim(ctr(m,1)+[-20, 20]);
%
%     subplot(3, 2, 3:4);
%     plot(neuron.C_raw(m, :));
%     xlim([0, 2000]);
%     subplot(3, 2, 5:6);
%     plot(neuron_ica.C_raw(match_cnmfe.ind_spatial(m), :));
%     disp([m, snr_cnmfe(m), snr_ica(match_cnmfe.ind_spatial(m))]);
%     xlim([0, 2000]);
%     pause;
% end
% ids_match = [5, 10, 13, 21, 24, 32, 38, 41, 89, 91, 105, 106, 290, 389];
% ids_match = [5, 9, 21, 25, 39, 54, 74, 106, 110, 146, 389];
% ids_match = linspace(1, K_match-50, 11);
ind = [39542, 50071, 26346, 8033, 24818, 59294, 31370, 33178, 16684, 28232, 30887];
[~, ids_ica] = max(neuron_ica.A(ind, :), [],2);
ids_match = match_ica.ind_max(ids_ica);
% ids_ica = match_cnmfe.ind_spatial(ids_match);
ctr = neuron.estCenter();
ctr_ica =  neuron_ica.estCenter();

%%
figure('papersize', [2, 5]);
init_fig;
for m=1:10
    axes('position', [0.505-0.5*mod(m,2), 1.003-ceil(m/2)*0.2, 0.49,0.196]);
    r0 = round(ctr(ids_match(m), 1));
    c0 = round(ctr(ids_match(m), 2));
    if r0<23
        r0 = 23;
    elseif r0>=d1-18
        r0 = d1-18;
    end
    if c0<21
        c0 = 21;
    elseif c0>=d2-20
        c0 = d2-20;
    end
    img = neuron.reshape(neuron.A(:, ids_match(m)),2);
    temp = img(r0+(-22:18), c0+(-20:20));
    imagesc(temp/max(temp(:)), [0,1]);
    
    axis equal off tight;
    text(1, 6 , num2str(m), 'color', 'w', 'fontweight', 'bold', 'fontsize', 20);
end
if export_fig
    saveas(gcf, sprintf('%s/match_spatial_cnmfe.fig', output_folder) );
    saveas(gcf, sprintf('%s/match_spatial_cnmfe.pdf', output_folder));
end

% PCA/ICA
figure('papersize', [2, 5]);
init_fig;
for m=1:10
    axes('position', [0.505-0.5*mod(m,2), 1.003-ceil(m/2)*0.2, 0.49,0.196]);
    r0 = round(ctr(ids_match(m), 1));
    c0 = round(ctr(ids_match(m), 2));
    if r0<23
        r0 = 23;
    elseif r0>=d1-18
        r0 = d1-18;
    end
    if c0<21
        c0 = 21;
    elseif c0>=d2-20
        c0 = d2-20;
    end
    img = neuron.reshape(A_ica_before_trim(:, ids_ica(m)), 2);
    temp = img(r0+(-22:18), c0+(-20:20));
    imagesc(temp/max(temp(:)), [-0.2,1]);
    
    axis equal off tight;
    text(1, 6, num2str(m), 'color', 'w', 'fontweight', 'bold', 'fontsize', 20);
end

if export_fig
    saveas(gcf, sprintf('%s/match_spatial_ica.fig', output_folder) );
    saveas(gcf, sprintf('%s/match_spatial_ica.pdf', output_folder));
end
%% temporal

figure('papersize', [6.5, 5]);
init_fig;
set(gcf, 'defaultAxesFontSize', 16);
col = colormap(cool);
C = neuron.C_raw(ids_match(1:10),:);
axes('position', [0.0, 0.0, 0.95, 1]); hold on;
for m=1:10
    y = C(m,:);
    plot((1:T)/neuron.Fs, y/max(y(1:1000))+10-m,'linewidth', 1, 'color', col(64-4*m,:));
    text(-8-1.5*(m==10), (10-m), num2str(m), 'fontsize', 20);
end
plot([88,97], [-0.5, -0.5], 'k', 'linewidth', 5);
% text(87, -1, '10 sec', 'fontsize', 18);
axis off;
axis([-10,100,-2, 10]);

if export_fig
    saveas(gcf, sprintf('%s/matched_temporal_cnmfe.fig', output_folder) );
    saveas(gcf, sprintf('%s/matched_temporal_cnmfe.pdf', output_folder));
end
% PCA/ICA
figure('papersize', [6.5, 5]);
init_fig;
set(gcf, 'defaultAxesFontSize', 16);
col = colormap(cool);
C = neuron_ica.C_raw(ids_match(1:10),:);
axes('position', [0.0, 0.0, 0.95, 1]); hold on;
for m=1:10
    y = C(m,:);
    plot((1:T)/neuron.Fs, y/max(y(1:1000))+10-m,'linewidth', 1, 'color', col(64-4*m,:));
    %     text(-8, (10-m), num2str(m), 'fontsize', 20);
end
plot([88,97], [-0.5, -0.5], 'k', 'linewidth', 5);
% text(87, -1, '10 sec', 'fontsize', 18);
axis off;
axis([-10,100,-2, 10]);
if export_fig
    saveas(gcf, sprintf('%s/matched_temporal_ica.fig', output_folder) );
    saveas(gcf, sprintf('%s/matched_temporal_ica.pdf', output_folder));
end
%% plot SNRs
figure('papersize', [3.5,5]);
init_fig;
set(gcf, 'defaultAxesFontSize', 14);
plot(snr_ica(1:K_match),snr_cnmfe(1:K_match),  '.k', 'markersize', 10);
hold on;
for m=1:length(ids_match)
    plot(snr_ica(ids_ica(m)), snr_cnmfe(ids_match(m)), 'o', ...
        'markersize', 10, 'markerfacecolor', col(64-4*m, :), 'markeredgecolor', 'none');
end
set(gca, 'position', [0.23, 0.15, 0.72, 0.83]);
hold on
axis  tight;

xlim([0, max(snr_ica(1:K_match))*1.1]);
ylim([0, max(snr_cnmfe(1:K_match))*1.1]);
ylabel('CNMF-E');
xlabel('PCA/ICA');
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
set(gca, 'xtick', 10.^(-1:1));
set(gca, 'ytick', 10.^(-1:2));

xlim([0.1, 10]);
ylim([0.1, 100]);
plot(get(gca, 'xlim'), get(gca, 'xlim'), '-r'); hold on;
if export_fig
    saveas(gcf, sprintf('%s/snr_pca_ica.fig', output_folder) );
    saveas(gcf, sprintf('%s/snr_pca_ica.pdf', output_folder));
end

%%
% show neuron contours
Coor = neuron.show_contours(0.6);

% create a video for displaying the
avi_filename = neuron.show_demixed_video(save_demixed, kt, [], amp_ac, range_ac, range_Y, multi_factor);
copyfile(avi_filename, demixed_video);

% save neurons shapes
neuron.save_neurons();

% %% save video
% kt = 4;     % play one frame in every kt frames
% save_avi = true;
% t_begin = 1;
% t_end = T;
%
% center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2))*1.25;
% range_res = [-1,1]*center_ac;
% range_ac = center_ac+range_res;
% multi_factor = 8;
% center_Y = min(Y(:)) + multi_factor*center_ac;
% range_Y = center_Y + range_res*multi_factor;
% avi_filename = sprintf('../../Videos/striatum_decompose.avi');
% cnmfe_save_video;
%
%
% %% demixing video
% file_nm = sprintf('../../Videos/blood_vessel_demix.avi');
% neuron.playAC(file_nm, [], [t_begin, t_end]);
%








