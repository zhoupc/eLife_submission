%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
addpath('./extra'); 
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_BG_subfigs'];
results_folder = sprintf('%sfig_BG%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_bg.mat'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% parameters to use
K = 50;    % number of neurons
K_bg = 25;  %background sources
d1 = 256;   % height
d2 = 256;   % weight
gSig = 3;   % gaussian width of cell body
d_min = 5;  % minimum distance between neurons
mu_bar = 0.01;  % mean spike counts per frame
T = 1000;   % number of frames
seed = 2;
tau_d_bar = 6;  % mean decay time constant, unit: frame
tau_r_bar = 1;  % mean rising time constant, unit: frame
minIEI = 20;    % minimum interval between events, unit: frame
minEventSize = 1; % minimum event size
sn = 5;
lambda = 20;     % scale the amplitude of cellular signals
rr = 5*gSig;
bg_neuron_ratio = 5;
ssub = 1;       % spatial downsampling factor
save default_parameters;
neuron = Sources2D('d1', d1, 'd2', d2, 'gSiz', 5*gSig);
neuron.Fs = 20;

%% simulate data 
sim_data;

%find the example pixel and determine its neighbors
x0 = 129;
y0 = 130;

rsub = (-rr):(rr);      % row subscript
csub = rsub;      % column subscript
[cind, rind] = meshgrid(csub, rsub);
R = sqrt(cind.^2+rind.^2);
neigh_kernel = (R>=rr) .* (R<rr+1);
[r_shift, c_shift] = find(neigh_kernel);
r_shift = r_shift - rr -1;
c_shift = c_shift - rr - 1;

r_neigh = y0 + r_shift;
c_neigh = x0 + c_shift;
ind_neigh = sub2ind([d1,d2], r_neigh, c_neigh);
ind_center = sub2ind([d1,d2], y0, x0);

%% Figure 2B, BG model 
figure('papersize', [5, 5]);
init_fig;
axes('position', [0.01, 0.01, 0.98, 0.98]);
rng(1);
neuron.A = A;
img = neuron.overlapA(17, 0);
imagesc(65536-img);
axis equal tight; hold on;
plot(x0, y0, 'og', 'markersize', 10, 'markerfacecolor', 'g');
plot(c_neigh, r_neigh, '.k', 'markersize', 30);
xlim(x0+1.5*[-rr, rr]);
ylim(y0+1.5*[-rr, rr]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);

if export_fig
    saveas(gcf, sprintf('%s/example_circle.fig', output_folder));
    saveas(gcf, sprintf('%s/example_circle.pdf', output_folder));
end

%% show one example frame
% Figure 1E 
figure('papersize', [5, 5]);
init_fig;
neuron.image(Y(:, 20));
colormap gray
axis equal tight;
hold on;
% plot(x0, y0, '+r', 'markersize', 10);
% plot(c_neigh, r_neigh, '.b', 'markersize', 10);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
if export_fig
    saveas(gcf, sprintf('%s/example_frame.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end

% Figure 2A 
set(gca,'position', [0.01, 0.01, 0.98, 0.98]);
colormap parula;
plot(x0, y0, 'og', 'markersize', 10, 'markerfacecolor', 'g');
plot(c_neigh, r_neigh, '.k', 'markersize', 10);
if export_fig
    saveas(gcf, sprintf('%s/example_frame_neuron.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_neuron.pdf', output_folder));
end
%% Figure 1E
figure('papersize', [5, 5]);
init_fig;
neuron.A = A;
rng(1);
img_neurons = neuron.overlapA([], 0.2);
neuron.image(imcomplement(img_neurons));
axis equal tight; hold on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
if export_fig
    saveas(gcf, sprintf('%s/example_neurons.fig', output_folder));
    saveas(gcf, sprintf('%s/example_neurons.pdf', output_folder));
end

% show localized background sources
figure('papersize', [5, 5]);
init_fig;
neuron.A = G(:, 3:end);
img_localBG = neuron.overlapA([], 0.3);
neuron.image(imcomplement(img_localBG));
axis equal tight; hold on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
if export_fig
    saveas(gcf, sprintf('%s/example_localBG.fig', output_folder));
    saveas(gcf, sprintf('%s/example_localBG.pdf', output_folder));
end

% show blood vessel
figure('papersize', [5, 5]);
init_fig;
img_vessel = neuron.reshape(G(:,2), 2);
neuron.image(imcomplement(img_vessel));
axis equal tight; hold on;
set(gca, 'xtick', []);
set(gca, 'ytick', []); colormap gray;
if export_fig
    saveas(gcf, sprintf('%s/example_bloodVessel.fig', output_folder));
    saveas(gcf, sprintf('%s/example_bloodVessel.pdf', output_folder));
end


% show blood vessel
figure('papersize', [5, 5]);
init_fig;
img_global = neuron.reshape(G(:,1), 2);
neuron.image((img_global));
axis equal tight; hold on;
set(gca, 'xtick', []);
set(gca, 'ytick', []); colormap gray;
if export_fig
    saveas(gcf, sprintf('%s/example_global.fig', output_folder));
    saveas(gcf, sprintf('%s/example_global.pdf', output_folder));
end
% show calcium traces and run regression
y_neigh = Y(ind_neigh, :);
y_center = Y(ind_center, :);
y_neigh = bsxfun(@minus, y_neigh, mean(y_neigh,2));
y_center_mean = mean(y_center);
y_center = y_center - y_center_mean;
w = (y_neigh*y_neigh')\(y_neigh*y_center');
y_center_fit = w' * y_neigh;

%% Figure 2C 
neuron.Fs = 10;
t = (1:T)/neuron.Fs;
figure('papersize', [10,3]);
init_fig;
axes('position', [0.01, 0.01, 0.98, 0.98]); hold on;
for m=1:10:length(ind_neigh)
    h1 = plot(t, y_neigh(m,:) + m*20, 'k', 'linewidth', 0.5);
end
%set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
ylim([-200, 2500]);
box on;%xlabel('Frame #');
axis off;
plot([t(end)-10, t(end)], -150*[1,1], 'k', 'linewidth', 8);
if export_fig
    saveas(gcf, sprintf('%s/example_neighbors.fig', output_folder)); %#ok<*UNRCH>
    saveas(gcf, sprintf('%s/example_neighbors.pdf', output_folder));
end

h2 = plot(t, y_center+(m+20)*20, 'g', 'linewidth', 1);
% legend([h2, h1],{'center pixel', 'neighboring pixels'}, 'orientation',...
%     'horizental', 'position', [0.35, 0.92, 0.6, 0.05]);
% axis tight; box off;
if export_fig
    saveas(gcf, sprintf('%s/example_neighbors_center.fig', output_folder)); %#ok<*UNRCH>
    saveas(gcf, sprintf('%s/example_neighbors_center.pdf', output_folder));
end

%% Figure 2D 
figure('papersize', [18, 4]);
init_fig; axes('position', [0.01, 0.01, 0.98, 0.98]); hold on;
set(gcf, 'defaultAxesfontsize', 8);
plot(t, y_center+y_center_mean+190, '-g', 'linewidth', 2); alpha(.5); hold on;
axis tight;
plot(t, y_center_fit+y_center_mean-80, 'b', 'linewidth', 3);
plot(t, Bf(ind_center,:)+Bc(ind_center,:)-80, '-.y', 'linewidth', 2);
set(gca, 'xticklabel', []);
box on;
plot(t, y_center-y_center_fit-100, 'r', 'linewidth', 3);
plot(t, Ysignal(ind_center,:)-100, '-.c', 'linewidth', 2);
legend('fluo. trace', 'estimated BG','true BG', 'BG-subtracted trace', ...
    'true neural signal',  'orientation', 'horizental', 'location', 'northwest');
plot([t(end)-10, t(end)], -200*[1,1], 'k', 'linewidth', 8);
ylim([-250, 800]);
xlim([t(1)-1, t(end)+1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'fontsize', 35);
if export_fig
    saveas(gcf, sprintf('%s/example_background.fig', output_folder));
    saveas(gcf, sprintf('%s/example_background.pdf', output_folder));
end

%% different lines for Figure 2D 
figure('papersize', [10,1.8]);
init_fig;
plot(y_center+y_center_mean, '-k', 'linewidth', 1.5); alpha(.5); hold on;
axis tight;
set(gca, 'xticklabel', []);
axis tight;
legend('observed fluorescence');
box on;

if export_fig
    saveas(gcf, sprintf('%s/yi.fig', output_folder));
    saveas(gcf, sprintf('%s/yi.pdf', output_folder));
end

figure('papersize', [10,1.8]);
init_fig;
plot(y_center_fit+y_center_mean, 'b', 'linewidth', 1.5); hold on;
% plot(Bf(ind_center,:)+Bc(ind_center,:), '-.r', 'linewidth', 1.5);
axis tight;
set(gca, 'xticklabel', []);
axis tight;
legend('estimated background');
box on;

if export_fig
    saveas(gcf, sprintf('%s/bi.fig', output_folder));
    saveas(gcf, sprintf('%s/bi.pdf', output_folder));
end
%
figure('papersize', [10, 2.2]);
init_fig; hold on;
plot(y_center-y_center_fit, 'g');
plot(Ysignal(ind_center,:), '-.r');
axis tight;
xlabel('Frame #');
legend('background subtracted trace', 'truth', 'orientation', 'horizental');
ylim([-30, 170]);
box on;
if export_fig
    saveas(gcf, sprintf('%s/example_cellular.fig', output_folder));
    saveas(gcf, sprintf('%s/example_cellular.pdf', output_folder));
end

%% use CNMF-E to estimate the background activity
try
    Ybg = results_data.Ybg;
    Ybg_cnmfe = results_data.Ybg_cnmfe;
    RSS_cnmfe = results_data.RSS_cnmfe;
    CC_cnmfe = results_data.CC_cnmfe;
catch
    Ybg = neuron.localBG(Y, 1, rr, [], [], 5);
    Ybg_cnmfe = bsxfun(@minus, Ybg, mean(Ybg,2));
    RSS_cnmfe = sum((Ybg_cnmfe(:)-Bf(:)).^2);
    CC_cnmfe = sum(Ybg_cnmfe.*Bf,2)./std(Ybg_cnmfe, 0, 2)./std(Bf,0,2)/T;
    save(results_file, 'RSS_cnmfe', 'CC_cnmfe', 'Ybg_cnmfe', 'Ybg', '-append');
end
%% Figure 2E and Figure 2F
figure('papersize', [5, 5]);
init_fig;
neuron.image(Y(:,20)-Ybg(:,20), [5, 80]);% colorbar;
axis equal off tight;
if export_fig
    saveas(gcf, sprintf('%s/example_bg_cnmfe_subtract.fig', output_folder));
    saveas(gcf, sprintf('%s/example_bg_cnmfe_subtract.pdf', output_folder));
end

figure('papersize', [5, 5]);
init_fig;
neuron.image(Ysignal(:,20), [5, 80]);
axis equal off tight;
if export_fig
    saveas(gcf, sprintf('%s/example_ac.fig', output_folder));
    saveas(gcf, sprintf('%s/example_ac.pdf', output_folder));
end

%% use rank -1 NMF to estimate the background
% Figure 2G 
[W,R] = nnmf(Y, 1);
figure('papersize', [5, 5]);
init_fig;
neuron.image(Y(:,20)-W*R(:,20), [5, 80]);% colorbar;
axis equal off tight;
if export_fig
    saveas(gcf, sprintf('%s/example_bg_cnmf_subtract.fig', output_folder));
    saveas(gcf, sprintf('%s/example_bg_cnmf_subtract.pdf', output_folder));
end

%% save video
% S2 Video 
avi_filename = [video_folder, 'background_comparison.avi'];
if ~exist(avi_filename, 'file')
    save_video;
end

%% use high rank NMF to estiamte the background  
try 
    RSS_nmf = results_data.RSS_nmf; 
    CC_nmf = results_data.CC_nmf; 
    nmf_rank = results_data.nmf_rank; 
catch
    sim_data;
    nmf_rank = 1:2:30;
    Nsim = 10;
    RSS_nmf = zeros(10, length(nmf_rank));
    CC_nmf = zeros(d1*d2, 10, length(nmf_rank));
    B = Bf + Bc;
    std_B = std(Bf, 0,2);
    % Y = Y-Ysignal;
    for m=1:length(nmf_rank)
        tic;
        for n=1:Nsim %#ok<*PFTUS>
            [W,R] = nnmf(Y, nmf_rank(m));
            Ybg_nmf = W*R;
            RSS_nmf(n, m) = sum((Ybg_nmf(:)-B(:)).^2);
            Ybg_nmf = bsxfun(@minus, Ybg_nmf, mean(Ybg_nmf,2));
            CC_nmf(:,n,m) = sum(Ybg_nmf.*Bf,2)./std(Ybg_nmf, 0, 2)./std_B/T;
        end
        disp([m, toc]);
    end
    
    save(results_file, 'RSS_nmf', 'CC_nmf', 'nmf_rank', '-append'); 
end

Nsim = size(RSS_nmf,1);

%% Figure 2H 
figure('papersize', [4,4]);
init_fig;
set(gcf, 'defaultaxesfontsize', 20);axes('position', [0.25, 0.25, 0.7, 0.7]);
CC = squeeze(mean(CC_nmf,1));
plot([0, nmf_rank(end)+.5], mean(CC_cnmfe)*ones(1,2), 'r'); hold on;
temp = repmat(nmf_rank, Nsim,1);
plot(temp(:), CC(:), 'ok');
% errorbar(nmf_rank, mean(RSS_nmf,1), std(RSS_nmf,0,1), 'g');
xlim([1, nmf_rank(end)+0.5]);
ylim([0.79, 1]);
legend('CNMF-E', 'high rank NMF', 'orientation', 'vertical', ...
    'location', 'southwest');
set(gca, 'ytick', [0.8, 0.9, 1]);
xlabel('NMF rank');
ylabel('Mean Corr.');

if export_fig
    saveas(gcf, sprintf('%s/example_corr_nmf.fig', output_folder));
    saveas(gcf, sprintf('%s/example_corr_nmf.pdf', output_folder));
end

% %% run sc-CNMF
% if ~exist('results_data.mat', 'file')
%     disk_sizs =2:15;
%     sim_data;
%     CC_sc = zeros(d1*d2, length(disk_sizs));
%     RSS_sc = zeros(length(disk_sizs), 1);
%     std_B = std(Bf, 0,2);
%     for m=1:length(disk_sizs)
%         tic;
%         kernel = strel('disk', disk_sizs(m));
%         Ysc = neuron.reshape(zeros(size(Y)), 2);
%         for t=1:T
%             Ysc(:, :, t) = imopen(neuron.reshape(Y(:, t), 2), kernel);
%         end
%         Ysc = neuron.reshape(Ysc, 1);
%         Ybg_sc = bsxfun(@minus, Ysc, mean(Ysc,2));
%         CC_sc(:,m) = sum(Ybg_sc.*Bf,2)./std(Ybg_sc, 0, 2)./std_B/T;
%         Yres = bsxfun(@minus, Ysc-Bf, Bc);
%         RSS_sc(m) = sum(Yres(:).^2);
%         disp([m, disk_sizs(m), toc]);
%     end
%     
%     save results_data CC_sc RSS_sc disk_sizs -append
% else
%     load results_data;
% end
% 
% 
% figure('papersize', [6,6]);
% init_fig;
% set(gcf, 'defaultaxesfontsize', 20);axes('position', [0.25, 0.25, 0.7, 0.7]);
% CC = mean(CC_sc,1);
% plot([0, disk_sizs(end)*2+2], mean(CC_cnmfe)*ones(1,2), 'r'); hold on;
% plot(disk_sizs*2+1, CC, 'ok');
% xlim([0, disk_sizs(end)*2+2]);
% % ylim([0.79, 1]);
% legend('CNMF-E', 'sc-CNMF', 'orientation', 'vertical', ...
%     'location', 'southwest');
% % set(gca, 'ytick', [0.8, 0.9, 1]);
% xlabel('kernel size');
% ylabel('Mean Corr.');
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_corr_scCNMF.fig', output_folder));
%     saveas(gcf, sprintf('%s/example_corr_scCNMF.pdf', output_folder));
% end
% 
% figure('papersize', [6,6]);
% init_fig;
% set(gcf, 'defaultaxesfontsize', 20);axes('position', [0.25, 0.25, 0.7, 0.7]);
% plot([0, disk_sizs(end)*2+2], RSS_cnmfe*ones(1,2), 'r'); hold on;
% plot(disk_sizs*2+1, RSS_sc, 'ok');
% xlim([0, disk_sizs(end)*2+2]);
% % ylim([0.79, 1]);
% legend('CNMF-E', 'sc-CNMF', 'orientation', 'vertical', ...
%     'location', 'northwest');
% % set(gca, 'ytick', [0.8, 0.9, 1]);
% xlabel('kernel size');
% ylabel('RSS');
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_RSS_scCNMF.fig', output_folder));
%     saveas(gcf, sprintf('%s/example_RSS_scCNMF.pdf', output_folder));
% end
% 
% %% run suite2P
% ops.TileFactor = 1;
% Lx = d2;
% Ly = d1;
% type = 'raisedcosyne';
% ops.diameter = neuron.options.gSiz;
% 
% if ~exist('results_data.mat', 'file')
%     neuropilRange=0.25:0.25:4;
%     sim_data;
%     CC_suite2p = zeros(d1*d2, length(neuropilRange));
%     RSS_suite2p = zeros(length(neuropilRange), 1);
%     ind_neuron = (sum(A,2)>0);
%     std_B = std(Bf, 0,2);
%     for m=1:length(neuropilRange)
%         tic;
%         ops.neuropilRange = neuropilRange(m);
%         Sbg = getNeuropilBasis(ops, Ly, Lx, type);
%         n_bg = (Sbg'*Sbg)\(Sbg'*Y);
%         
%         Ybg_suite2p = Sbg*bsxfun(@minus, n_bg, mean(n_bg,2));
%         CC_suite2p(:,m) = sum(Ybg_suite2p.*Bf,2)./std(Ybg_suite2p, 0, 2)./std_B/T;
%         Yres = bsxfun(@minus, Sbg*n_bg-Bf, Bc);
%         RSS_suite2p(m) = sum(Yres(:).^2);
%         ResStd_suite2p(:,m) = std(Yres, 0, 2);
%         ResMean_suite2p(:, m) = mean(Yres,2);
%         disp([m, neuropilRange(m), toc]);
%     end
%     
%     save results_data CC_suite2p RSS_suite2p ResMean_suite2p ResStd_suite2p ind_neuron neuropilRange -append
% else
%     load results_data;
% end
% 
% 
% figure('papersize', [6,6]);
% init_fig;
% set(gcf, 'defaultaxesfontsize', 20);axes('position', [0.25, 0.25, 0.7, 0.7]);
% CC = mean(CC_suite2p,1);
% plot([0, neuropilRange(end)+.5], mean(CC_cnmfe)*ones(1,2), 'r'); hold on;
% plot(neuropilRange, CC, 'ok');
% xlim([0, neuropilRange(end)+0.5]);
% % ylim([0.79, 1]);
% legend('CNMF-E', 'Suite2P', 'orientation', 'vertical', ...
%     'location', 'southwest');
% % set(gca, 'ytick', [0.8, 0.9, 1]);
% xlabel('NeuropilRange');
% ylabel('Mean Corr.');
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_corr_suite2p.fig', output_folder));
%     saveas(gcf, sprintf('%s/example_corr_suite2p.pdf', output_folder));
% end
% 
% 
% figure('papersize', [6,6]);
% init_fig;
% set(gcf, 'defaultaxesfontsize', 20);axes('position', [0.25, 0.25, 0.7, 0.7]);
% plot([0, neuropilRange(end)+.5], RSS_cnmfe*ones(1,2), 'r'); hold on;
% plot(neuropilRange, RSS_suite2p, 'ok');
% xlim([0, neuropilRange(end)+0.5]);
% % ylim([0.79, 1]);
% legend('CNMF-E', 'Suite2P', 'orientation', 'vertical', ...
%     'location', 'northwest');
% set(gca, 'ytick', [0.8, 0.9, 1]);
% xlabel('NeuropilRange / Neuron diameter');
% ylabel('Mean Corr.');
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_rss_suite2p.fig', output_folder));
%     saveas(gcf, sprintf('%s/example_rss_suite2p.pdf', output_folder));
% end

%%
%% performances of the ring model and NMF model on different number of the background sources 
try
    CC_all = results_data.CC_all;
    CC_all_nmf = results_data.CC_all_nmf;
    RSS_all = results_data.RSS_all;
    RSS_all_nmf = results_data.RSS_all_nmf;
    K_bg_all = results_data.K_bg_all;
catch
    load default_parameters;
    neuron = Sources2D('d1', d1, 'd2', d2, 'gSiz', 5*gSig);
    sim_data;
    
    Nval = 10;      % number of possible values for the background rank
    K_bg_all = unique(ceil(logspace(0,2,Nval)));
    Nval = length(K_bg_all);
    RSS_all = zeros(Nval,1);
    CC_all = zeros(d1*d2,Nval,1);
    
    RSS_all = zeros(Nval,1);
    CC_all = zeros(d1*d2,Nval,1);
    
    RSS_all_thresh = RSS_all;
    CC_all_thresh = CC_all;
    
    parfor m=1:Nval
        tic;
        K_bg = K_bg_all(m);
        % generate new background
        [G, G_coor] = gen_spatial(K_bg, [d1,d2], gSig*5, d_min, seed*2);
        % temporal
        F = zeros(K_bg_all(m), T);
        for n=1:K_bg_all(m)
            g = rand(1)*0.3+0.5;
            tmp_f = randn(1, T);
            for t=2:T
                tmp_f(t) = tmp_f(t-1)*g + tmp_f(t);
            end
            tmp_f = tmp_f - min(tmp_f);
            F(n,:) = tmp_f/std(tmp_f);
        end
        
        % replace the first background source with a global fluctuations
        G(:, :,1) = imresize(G_global, [d1, d2]);
        F(1,:) = resample(F_global, T, length(F_global));
        
        % replace the second background source with a blood vessels
        if K_bg>1
            img = zeros(d1,d2);
            x = (1:d2);
            y =  d1 - round(4*(x/d2-1/2).^3*d1 + d1/2);
            y(y<1) = 1; y(y>d1) = d1;
            ind = sub2ind([d1,d2], y, x);
            img(ind) = 1;
            img = imfilter(img, fspecial('gaussian', 10, gSig));
            G(:, :, 2) = img/max(img(:));
            F(2,:) = resample(F_vessel, T, length(F_vessel));
        end
        
        % change the amplitude of each background source
        ind_center = sub2ind([d1,d2], round(G_coor.y0), round(G_coor.x0));
        temp = G(:,:,1);
        if K_bg == 1
            F = F*50;
        elseif K_bg == 2
            F = bsxfun(@times, F, 10*[5;2]);
        else
            F = bsxfun(@times, F, 10*[5; 2; 3* sqrt(temp(ind_center(3:end)))]);
        end
        G = reshape(G, d1*d2,[]);
        
        % combine results and generate the observed fluorescence
        Fmean = mean(F,2);
        tmp_Bf = G * bsxfun(@minus, F, Fmean);
        tmp_Bc = G * Fmean * ones(1,T);
        Y = tmp_Bf+tmp_Bc+Ysignal+sn*E;
        std_Bf = std(tmp_Bf,0,2);
        
        %         run CNMF-E
        Ybg = local_background(reshape(Y, d1, d2, T), ssub, rr);
        Ybg = reshape(Ybg, [], T);
        RSS_all(m) = sum((Ybg(:)-tmp_Bf(:)-tmp_Bc(:)).^2);
        Ybg = bsxfun(@minus, Ybg, mean(Ybg,2));
        CC_all(:,m) = mean(Ybg.*tmp_Bf,2)./std(Ybg, 0, 2)./std_Bf;
        
        %         run rank-1 NMF
        [W,R] = nnmf(Y, 1);
        Ybg = W*R;
        RSS_all_nmf(m) = sum((Ybg(:)-tmp_Bf(:)-tmp_Bc(:)).^2);
        Ybg = bsxfun(@minus, Ybg, mean(Ybg,2));
        CC_all_nmf(:,m) = mean(Ybg.*tmp_Bf,2)./std(Ybg, 0, 2)./std_Bf;
        
        fprintf('%d: background rank %d, time ellapsed %.3f \n', m, K_bg, toc);
    end
    save(results_file, 'CC_all', 'CC_all_nmf', 'RSS_all', 'RSS_all_nmf',  'K_bg_all', '-append') ;
end

%% plot results. Figure 2I 
figure('papersize', [4,4]);
init_fig;
set(gcf, 'defaultaxesfontsize', 20);
axes('position', [0.23, 0.25, 0.7, 0.7]);
plot(K_bg_all, mean(CC_all,1), '-or');
hold on;
plot(K_bg_all,mean(CC_all_nmf,1), '-dk');
set(gca, 'xscale', 'log');
set(gca, 'xtick', [2, 10, 100]);
xlim([2, K_bg_all(end)+0.5]);
ylim([0.79, 1.0]);
set(gca, 'ytick', [0.8, 0.9, 1.0]);
legend('CNMF-E', 'rank-1 NMF', 'orientation', 'vertical', 'location', ...
    'southwest');
% title('mean correlation');
xlabel('# of BG sources');
ylabel('Mean corr.');
% if export_fig
%     saveas(gcf, sprintf('%s/example_RSS_corr_high_rank.fig', output_folder));
%     saveas(gcf, sprintf('%s/example_RSS_corr_high_rank.pdf', output_folder));
% end
if export_fig
    saveas(gcf, sprintf('%s/example_corr_high_rank.fig', output_folder));
    saveas(gcf, sprintf('%s/example_corr_high_rank.pdf', output_folder));
end

%% remove the added path 
rmpath('../functions/');
rmpath('./extra'); 