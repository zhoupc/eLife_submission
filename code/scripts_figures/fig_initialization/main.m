%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_INIT_subfigs'];
results_folder = sprintf('%sfig_initialization%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_initialization.mat'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% load the background and noise levels extracted from the a typical microendoscopic data
load data_BG;
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
neuron.updateParams('gSig', gSig, 'gSiz', gSiz);

%% simulate cellular signals
sim_AC;
ind_center = sub2ind([d1,d2], round(coor.y0), round(coor.x0));
A = bsxfun(@times, A, reshape(sn(ind_center), 1, []));   %adjust the absolute amplitude neuron signals based on its spatial location
%% simulate white noise
E = bsxfun(@times, randn(d1*d2, T), sn);

%% generate the observed data
Bc = D*mean(F,2);       % time-invariant baseline
Bf = D*bsxfun(@minus, F, mean(F,2));  % fluctuating background
Y = neuron_amp*A*C + bsxfun(@plus, Bf, Bc) + E;  % observed fluorescence

%% example frame. Figure 3A 
% pick an isolated neuron as example
cell_id1 = 99;
center_id1 = ind_center(cell_id1);
cell_id2 = [139, 152];
center_id2 = ind_center(cell_id2);

figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
ind_frame = 578;
neuron.image(Y(:, ind_frame), [400, 2200]);
set(gca,'position',[0 0 1 1],'units','normalized')
colormap gray;
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

if export_fig
    saveas(gcf, sprintf('%s/example_frame.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end
%% example of single neuron. Figure 3D
% plot spatial
figure('papersize', [3,3]);
init_fig;
img = zeros(d1, d2, 3);
img(:,:,3) = neuron.reshape(A(:, cell_id1), 2);
img(img==0) = min(img(img~=0)/2);
img = uint16(img/max(img(:))*65536);
imagesc(65536-img);  hold on;
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot(x0, y0, '*b', 'markersize', 18);
ind_0 = sub2ind([d1,d2], y0, x0);
x1 = x0+gSig;
y1 = y0+gSig;
plot(x1, y1, '+m', 'markersize', 18);
ind_1 = sub2ind([d1,d2], y1, x1);
x2 = x0+2*gSig;
y2 = y0+2*gSig;
plot(x2, y2, 'xc', 'markersize', 18);
ind_2 = sub2ind([d1,d2], y2, x2);
axis equal tight; box on; 
set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
xlim(x0+[-1,1]*gSiz);
ylim(y0+[-1,1]*gSiz);

if export_fig
    saveas(gcf, sprintf('%s/single_neuron_spatial.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/single_neuron_spatial.pdf', output_folder));
end

% plot traces
figure('papersize', [8, 2.8]);
init_fig;

axes('position', [0.01, 0.0, 0.98, 1]); cla; 
y_0 = Y(ind_0,:); 
y_1 = Y(ind_1,:); 
y_2 = Y(ind_2,:); 

plot(y_0-mean(y_0)+600, 'b'); hold on; 
plot(-100, 600, '*b', 'markersize', 25); 
plot(-100, 300, '+m', 'markersize', 25); 
plot(-100, 0, 'xc', 'markersize', 25); 
plot(y_1-mean(y_1)+300, 'm'); 
plot(y_2-mean(y_2), 'c'); 
axis  off; 
xlim([-120, length(y_0)+5]); 
ylim([-200, 1500]); 
if export_fig
    saveas(gcf, sprintf('%s/single_neuron_temporal_raw.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/single_neuron_temporal_raw.pdf', output_folder));
end

% plot filtered results
% create spatial filter
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;
[nr, nc] = size(psf);
[indc, indr] = meshgrid((1:nc)-ceil(nc/2),(1:nr)-ceil(nr/2));
ind_shift = d1*(indc-1)+indr; % location of the surrounding pixels relative to the center pixel

figure('papersize', [8, 2.8]);
init_fig;

pos = [0.01, 0.0, 0.98, 1];
axes('position',pos);
temp = psf(:)'*Y(ind_0+ind_shift(:),:);
temp_baseline = median(temp);
true_c = C(cell_id1,:);
temp = psf(:)'*Y(ind_0+ind_shift(:),:)+6; 
plot(temp, 'b');
hold on; 
x = polyfit(true_c, temp, 1); 
plot(x(1)*true_c+x(2), '-.y');
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;
% text(-200, 10, '*', 'color', 'b', 'fontsize', 50);

plot(psf(:)'*Y(ind_1+ind_shift(:),:)+3, 'm');
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;

plot(psf(:)'*Y(ind_2+ind_shift(:),:), 'c');
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;

spk = find(diff(true_c)>1); 
for m=1:length(spk); 
    plot(spk(m), 1, 'vr');
end 
xlim([-120, length(y_0)+5]); 
ylim([-5, 25]); 
if export_fig
    saveas(gcf, sprintf('%s/single_neuron_temporal_filtered.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/single_neuron_temporal_filtered.pdf', output_folder));
end
%% two close neurons. Figure 3E 
figure('papersize', [3,3]);
init_fig;
img = zeros(d1, d2, 3);
img(:, :, 1:2) = neuron.reshape(A(:, fliplr(cell_id2)), 2);
img(:, :, 3) = sum(img, 3); 
img = uint16(img/max(img(:))*65536);
imagesc(65536-img);
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
ind_0 = sub2ind([d1,d2], y0, x0);
ind_1 = sub2ind([d1,d2], coor.y0(cell_id2(1)), coor.x0(cell_id2(1)));
ind_2 = sub2ind([d1,d2], coor.y0(cell_id2(2)), coor.x0(cell_id2(2)));
axis equal  tight;  hold on;
xlim(x0+[-1,1]*gSiz);
ylim(y0+[-1,1]*gSiz);
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
box on; set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
plot(x0, y0, '*b', 'markersize', 18);
plot(coor.x0(cell_id2(2)), coor.y0(cell_id2(2)), '+m', 'markersize', 18);
plot(coor.x0(cell_id2(1)), coor.y0(cell_id2(1)), 'xc', 'markersize', 18);

if export_fig
    saveas(gcf, sprintf('%s/double_neuron_spatial.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/double_neuron_spatial.pdf', output_folder));
end

% plot traces
figure('papersize', [8, 2.8]);
init_fig;

axes('position', [0.01, 0.0, 0.98, 1]); cla; 
y_0 = Y(ind_0,:); 
y_1 = Y(ind_1,:); 
y_2 = Y(ind_2,:); 

dh = 600; 
plot(y_0-mean(y_0)+2*dh, 'b'); hold on; 
plot(-100, 2*dh, '*b', 'markersize', 25); 
plot(-100, dh, '+m', 'markersize', 25); 
plot(-100, 0, 'xc', 'markersize', 25); 
plot(y_1-mean(y_1)+dh, 'm'); 
plot(y_2-mean(y_2), 'c'); 
axis  off; 
xlim([-120, length(y_0)+5]); 
ylim([-400, 2500]); 

if export_fig
    saveas(gcf, sprintf('%s/double_neuron_temporal_raw.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/double_neuron_temporal_raw.pdf', output_folder));
end

% filtered traces 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;
[nr, nc] = size(psf);
[indc, indr] = meshgrid((1:nc)-ceil(nc/2),(1:nr)-ceil(nr/2));
ind_shift = d1*(indc-1)+indr; % location of the surrounding pixels relative to the center pixel

figure('papersize', [8, 2.8]);
init_fig;
dh = 15; 
pos = [0.01, 0.0, 0.98, 1];
axes('position',pos);

temp = psf(:)'*Y(ind_0+ind_shift(:),:);
temp_baseline = median(temp);
true_c = C(cell_id2,:);
temp = psf(:)'*Y(ind_0+ind_shift(:),:)+2*dh; 
plot(temp, 'b'); 
hold on;
% plot(true_c/sum(true_c) * sum(temp-temp_baseline)+temp_baseline+10, 'r', 'linewidth', 1);
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;
% text(-200, 10, '*', 'color', 'b', 'fontsize', 50);

temp = psf(:)'*Y(ind_1+ind_shift(:),:)+dh; 
x = polyfit(true_c(1,:), temp, 1); 
plot(temp, 'm'); hold on; 
plot(x(1)*true_c(1,:)+x(2), '-.g'); 
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;

temp = psf(:)'*Y(ind_2+ind_shift(:),:); 
plot(temp, 'c'); 
x = polyfit(true_c(2,:)', temp', 1); 
plot(x(1)*true_c(2,:)+x(2), '-.r'); 
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
axis  off;

spk1 = find(S(cell_id2(1),:)); 
[~,ind_max]=max(S(cell_id2(1),:)); 
for m=1:length(spk1); 
    plot(spk1(m), 2*dh, '^m');
end 
spk2 = find(S(cell_id2(2),:)); 
for m=1:length(spk2); 
    plot(spk2(m), 2*dh, '^c');
end 
plot(ind_max, 5, 'vm'); 
xlim([-120, length(y_0)+5]); 
ylim([-5, 55]); 
if export_fig
    saveas(gcf, sprintf('%s/double_neuron_temporal_filtered.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/double_neuron_temporal_filtered.pdf', output_folder));
end

%% correlation image of the raw data 
ymean = mean(Y, 1);
ymean = ymean - mean(ymean);
ymean = ymean / norm(ymean,2);
temp = Y - Y*(ymean')*ymean;
Cn_raw = correlation_image(temp, 8, d1, d2);
% Cn_raw = correlation_image(Y, 8, d1, d2);  % subtract the largest
% component
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Cn_raw, [0, 1]);
set(gca,'position',[0 0 1 1],'units','normalized')
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

if export_fig
    saveas(gcf, sprintf('%s/example_Cn_raw.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_Cn_raw.pdf', output_folder));
end
%% compute correlation image and PNR
try 
    Cn = results_data.Cn; 
    PNR = results_data.PNR; 
catch
    [Cn, PNR] = neuron.correlation_pnr(Y);
    save(results_file, 'Cn', 'PNR', '-append');
end
% correlation image. Figure 3C
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Cn, [0, 1]);
set(gca,'position',[0 0 1 1],'units','normalized')
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

if export_fig
    saveas(gcf, sprintf('%s/example_Cn.fig', output_folder), 'psc2'); %#ok<*UNRCH>
    saveas(gcf, sprintf('%s/example_Cn.pdf', output_folder));
end
%%  peak to noise ratio
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
neuron.image(PNR);
set(gca,'position',[0 0 1 1],'units','normalized')
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

saveas(gcf, sprintf('%s/example_pnr.fig', output_folder), 'psc2');
saveas(gcf, sprintf('%s/example_pnr.pdf', output_folder));

%% Cn.*pnr
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Cn.*PNR);
set(gca,'position',[0 0 1 1],'units','normalized')
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

saveas(gcf, sprintf('%s/example_cn_pnr.fig', output_folder), 'psc2');
saveas(gcf, sprintf('%s/example_cn_pnr.pdf', output_folder));

%% average neuron activity. Figure 3B 
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
ind_frame = 578;
neuron.image(A*mean(C,2));
set(gca,'position',[0 0 1 1],'units','normalized')
axis equal off tight;
hold on;
% draw a square surrounding the first neuron
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g');
% draw a square surrounding the second neuron
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'r');

saveas(gcf, sprintf('%s/example_mean.fig', output_folder), 'psc2');
saveas(gcf, sprintf('%s/example_mean.pdf', output_folder));

%% run initialization
rng(1);
debug_on = false;
save_avi = [video_folder, 'sim_initialization.avi'];
patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = 0.9;     % minimum local correlation for a seeding pixel
min_pnr = 30;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd, 'gSig', gSig, 'gSiz', gSiz, 'deconv_flag', true);
neuron.options.deconv_options = struct('type', 'ar2', ...
    'method', 'thresholded', ...
    'optimize_pars', false, ...
    'optimize_smin', false);
% greedy method for initialization
tic;
[center, Cn, ~] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

neuron_init = neuron.copy;
save(results_file, 'neuron_init', 'D', 'F', 'A', 'C', 'sn', '-append'); 

%% show results
figure('papersize', [d2, d1]/max(d1,d2)*5);
init_fig;
set(gca, 'position', [0, 0, 1, 1]);
if ~exist('init_contour', 'var')
    init_contour = plot_contours(neuron.A, Cn, 0.95, 0, [], [], 2, 'k');
else
    plot_contours(neuron.A, Cn, 0.95, 0, [], init_contour, 2, 'k');
end
hold on;
colormap gray; 
plot(coor.x0, coor.y0, '.b', 'markersize', 10);
% plot(init_contour{end}(1,3:end), init_contour{end}(2,3:end), 'm', 'linewidth', 2);

if export_fig
    saveas(gcf, sprintf('%s/example_init.fig', output_folder));
    saveas(gcf, sprintf('%s/example_init.pdf', output_folder));
end
%% spatial components of three selected neurons
[~, ind] = max(neuron.A'*A(:,[cell_id1, cell_id2]), [], 1);   % find the best match

figure('papersize', [3,3]);
init_fig;
img = zeros(d1, d2, 3);
img(:,:,3) = neuron.reshape(neuron.A(:,ind(1)), 2);
img(img==0) = min(img(img~=0)/2);
img = uint16(img/max(img(:))*65536);
imagesc(65536-img);  hold on;
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
x0 = round(coor.x0(cell_id1));
y0 = round(coor.y0(cell_id1));
% plot(x0, y0, '*b', 'markersize', 18);
ind_0 = sub2ind([d1,d2], y0, x0);
x1 = x0+gSig;
y1 = y0+gSig;
% plot(x1, y1, '+m', 'markersize', 18);
ind_1 = sub2ind([d1,d2], y1, x1);
x2 = x0+2*gSig;
y2 = y0+2*gSig;
% plot(x2, y2, 'xc', 'markersize', 18);
ind_2 = sub2ind([d1,d2], y2, x2);
axis equal  tight;
xlim(x0+[-1,1]*gSiz);
ylim(y0+[-1,1]*gSiz);
box on; 
set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
if export_fig
    saveas(gcf, sprintf('%s/single_neuron_spatial_initialized.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/single_neuron_spatial_initialized.pdf', output_folder));
end
%% two close neurons
figure('papersize', [3,3]);
init_fig;
img = zeros(d1, d2, 3);
img(:, :, 1:2) = neuron.reshape(neuron.A(:, fliplr(ind(2:3))), 2);
img(:, :, 3) = sum(img,3); 
img = uint16(img/max(img(:))*65536);
imagesc(65536-img);
x0 = round(mean(coor.x0(cell_id2)));
y0 = round(mean(coor.y0(cell_id2)));
ind_0 = sub2ind([d1,d2], y0, x0);
ind_1 = sub2ind([d1,d2], coor.y0(cell_id2(1)), coor.x0(cell_id2(1)));
ind_2 = sub2ind([d1,d2], coor.y0(cell_id2(2)), coor.x0(cell_id2(2)));
axis equal  tight;  hold on;
xlim(x0+[-1,1]*gSiz);
ylim(y0+[-1,1]*gSiz);
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
box on; 
set(gca, 'xtick', []); 
set(gca, 'ytick', []); 

% % plot(x0, y0, '*b', 'markersize', 18);
% plot(coor.x0(cell_id2(2)), coor.y0(cell_id2(2)), 'xm', 'markersize', 18);
% plot(coor.x0(cell_id2(1)), coor.y0(cell_id2(1)), '+c', 'markersize', 18);

if export_fig
    saveas(gcf, sprintf('%s/double_neuron_spatial_initialized.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/double_neuron_spatial_initialized.pdf', output_folder));
end
%% spatial and temporal similarities  between the initialization and the ground truth
% match CNMF-E initialization
match_cnmfe_init = pair_neurons(A, C, neuron_init.A, neuron_init.C_raw);
ind_cnmfe_init = (match_cnmfe_init.ind_spatial == match_cnmfe_init.ind_temporal);

figure('papersize', [7, 6]);
init_fig;
set(gcf, 'defaultaxesfontsize', 30); 
axes('position', [0.23, 0.2, 0.75, 0.75]); 
plot(match_cnmfe_init.max_spatial(ind_cnmfe_init), match_cnmfe_init.max_temporal(ind_cnmfe_init), '+b');
hold on;
%plot(max_A, C_similar_raw(ind_max), '*g');
% set(gca, 'xtick', [0, 0.5, 1]);
% set(gca, 'ytick', [0, 0.5, 1]);
xlim([0.88, 1]);
ylim([0.88, 1]);
set(gca, 'ytick', 0.9:0.05:1); 
set(gca, 'xtick', 0.8:0.05:1); 
xlabel('Spatial similarity');
ylabel('Temporal similarity');
if export_fig
    saveas(gcf, sprintf('%s/example_init_similarity.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_init_similarity.pdf', output_folder));
end
%% check the sensitivity to the selection of min_pnr and min_corr
Cn_max = ordfilt2(Cn, gSig^2, true(gSig));
PNR_max = ordfilt2(PNR, gSig^2, true(gSig));
ind_max = or(Cn==Cn_max, PNR==PNR_max);
ind_true = neuron.reshape(sum(bsxfun(@gt, A, max(A,[],1)*exp(-1)), 2)>0, 2);
figure('papersize', [7, 6]);
init_fig;
axes('position', [0.23, 0.2, 0.75, 0.75]); 
% plot dots corresponding to neurons
plot(Cn(and(ind_max, ind_true)), PNR(and(ind_max, ind_true)), 'ob');
hold on;
plot(Cn(and(ind_max, ~ind_true)), PNR(and(ind_max, ~ind_true)), '*g');
plot(get(gca, 'xlim'), [1,1]*min_pnr, 'r');
plot([1,1]*min_corr, get(gca, 'ylim'), 'r');
xlim([0, 1]); 
ylim([0, 100]); 
x = legend('Pixels in neuron centers', 'Others', 'Thresholds', 'location', ...
    'northwest');
set(x, 'fontsize', 25);
xlabel('Local correlation');
ylabel('Peak-to-noise ratio');
set(gca, 'fontsize', 30); 
% plot dots corresponding to noises

if export_fig
    saveas(gcf, sprintf('%s/example_min_cn_pnr.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_min_cn_pnr.pdf', output_folder));
end

