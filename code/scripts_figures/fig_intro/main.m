%% initialize working space and prepare for the computation
clear; clc; close all;
addpath('../functions/');
work_dir = fileparts(mfilename('fullpath'));
prepare_env;
output_folder = [fig_folder, filesep, 'Fig_INTRO_subfigs'];
results_folder = sprintf('%sfig_intro%s', results_folder, filesep);
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end
results_file = [results_folder, 'fig_intro.mat'];

if ~exist(results_file, 'file')
    save(results_file, 'hash_cnmfe', '-v7.3');
end
results_data = matfile(results_file);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
export_fig = true;

%% load data 
nam = [data_folder, 'blood_vessel_10Hz.mat']; 
data = matfile(nam); 

Y = data.Y(:, :, 1:3000); 
[d1,d2,T] = size(Y);    % video dimension 
Y = double(Y);
neuron = Sources2D('d1', d1, 'd2', d2); 
Y = neuron.reshape(Y,1); 
Ymean = mean(Y,2); 
Ymedian = median(Y,2); 
neuron.Fs = 10; 

%% save an example video. S1 Video 
avi_name = [video_folder, 'example_microendoscopic_data.avi'];
if ~exist('avi_name', 'file')
    neuron.playMovie(Y, [1100, 2800], 'gray', avi_name, 1/neuron.Fs);
end

%% choose one example frame. Figure 1A 
ind_frame = 1000; 
x0 = 118;   % center of the selected neuron 
y0 = 56; 
gSiz = 11; 
dmin = round(gSiz/2); 

theta = linspace(0, 2*pi, 100); 
x1 = x0+cos(theta)*dmin; 
y1 = y0+sin(theta)*dmin; 

x2 = x0+[cos(theta)*dmin*1.2, cos(theta)*dmin*1.8];
y2 = y0+[sin(theta)*dmin*1.2, sin(theta)*dmin*1.8];

% raw image 
figure('papersize', [6, 5*d1/d2]);
init_fig;
axes('position', [0, 0, 0.9, 1]); 
neuron.image(Y(:, ind_frame));
colormap gray
axis equal tight;
hold on;
% plot(x0, y0, '+r', 'markersize', 10);
% plot(c_neigh, r_neigh, '.b', 'markersize', 10);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
    [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'g', 'linewidth', 4);
% fill(x1, y1, 'r', 'facealpha', 0.3, 'edgecolor', 'none'); 
% fill(x2, y2, 'g', 'facealpha', 0.3, 'edgecolor', 'none'); 

colorbar; 
if export_fig
    saveas(gcf, sprintf('%s/example_frame.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end

% raw image figure('papersize', [5, 5*d1/d2]);

figure('papersize', [6, 5*d1/d2]);
init_fig;
axes('position', [0, 0, 0.9, 1]); 
neuron.image(Y(:, ind_frame)-Ymean);
colormap gray
axis equal tight;
hold on;
% plot(x0, y0, '+r', 'markersize', 10);
% plot(c_neigh, r_neigh, '.b', 'markersize', 10);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
fill(x1, y1, 'm', 'facealpha', 0.3, 'edgecolor', 'none'); 
fill(x2, y2,'b', 'facealpha', 0.3, 'edgecolor', 'none'); 

% plot([x0-gSiz, x0+gSiz, x0+gSiz, x0-gSiz, x0-gSiz], ...
%     [y0-gSiz, y0-gSiz, y0+gSiz, y0+gSiz, y0-gSiz], 'b', 'linewidth', 4);
colorbar; 
if export_fig
    saveas(gcf, sprintf('%s/example_frame_no_baseline.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_no_baseline.pdf', output_folder));
end

%% compute correlation image. Figure 1B 
Cn = correlation_image(Y, 8, d1, d2); 
figure('papersize', [5, 4*d1/d2]); 
init_fig; 
axes('position', [0, 0, 0.9,1]); 
imagesc(Cn, [0.8, 1]); colorbar;  
axis equal off tight; 
if export_fig
    saveas(gcf, sprintf('%s/correlation_image.fig', output_folder));
    saveas(gcf, sprintf('%s/correlation_image.pdf', output_folder));
end

%% crop a small region 
Y = neuron.reshape(Y, 2); 
Ybox = Y(y0+(-gSiz:gSiz), x0+(-gSiz:gSiz),:); 
Ybox_mean = mean(Ybox, 3); 
Ybox_median = median(Ybox, 3); 

% raw image 
figure('papersize', [4,4]);
init_fig;
axes('position', [0, 0, 1, 1]); 
imagesc(Ybox(:,:,ind_frame)); 
colormap gray
axis equal tight off;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_box.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_box.pdf', output_folder));
end
% raw image without baseline 
figure('papersize', [4,4]);
init_fig;
axes('position', [0, 0, 1, 1]); 
imagesc(Ybox(:,:,ind_frame)-Ybox_median); 
colormap gray
axis equal tight off;
hold on;
if export_fig
    saveas(gcf, sprintf('%s/example_frame_box_no_baseline.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_box_no_baseline.pdf', output_folder));
end

%% Figure 1C
theta = linspace(0, 2*pi, 100); 
x1 = gSiz+1+cos(theta)*dmin; 
y1 = gSiz+1+sin(theta)*dmin; 
figure('papersize', [4,4]); 
init_fig; 

img = uint16(zeros(2*gSiz+1, 2*gSiz+1, 3));
Ybox = reshape(Ybox, 2*gSiz+1, 2*gSiz+1, T);
temp = Ybox(:,:, ind_frame)- Ybox_mean;

axes('position', [0,0,1,1]); 
imagesc(temp); colormap gray; 
hold on; 
fill(x1, y1, 'm', 'facealpha', 0.5, 'edgecolor', 'none'); 

x1 = gSiz+1+[cos(theta)*dmin*1.2, cos(theta)*dmin*1.8];
y1 = gSiz+1+[sin(theta)*dmin*1.2, sin(theta)*dmin*1.8];
fill(x1, y1, 'b', 'facealpha', 0.5, 'edgecolor', 'none'); 

if export_fig
    saveas(gcf, sprintf('%s/example_box_roi.fig', output_folder));
    saveas(gcf, sprintf('%s/example_box_roi.pdf', output_folder));
end

%% estimate the calcium traces. Figure 1D 
 Y = neuron.reshape(Y,1); 
[x, y] = meshgrid(-gSiz:gSiz, -gSiz:gSiz); 
dist_pixel = sqrt(x.^2+ y.^2); 
ind_center =  (dist_pixel<dmin); 
ind_boundary = (dist_pixel>dmin*1.2) & (dist_pixel<1.8*dmin); 
Ybox = reshape(Ybox, [], T); 
y1 = mean(Ybox(ind_center,:), 1)-mean(Ybox_median(ind_center));
y2 = mean(Ybox(ind_boundary(:),:)) -mean(Ybox_median(ind_boundary));
ci = y1-corr(y1', y2')*y2; 
t = (1:T)/neuron.Fs;

figure('papersize', [12.4, 3]); 
init_fig; 
set(gcf, 'defaultAxesFontSize', 20); 
axes('position', [0.01, 0.01, 0.98, 0.98]); 
dh = 250; 
plot(t, y1+2*dh, 'm'); hold on; 
plot(t, y2+dh, 'b'); 
plot(t, y1-y2, 'r'); 
legend('center ROI', 'periphery ROI', 'difference', 'orientation', 'horizental'); 
plot([260, 300], [-60,-60], 'k', 'linewidth', 6); 
axis([0, 300, -100, 3.5*dh]); 
plot([-5, -5, 305, 305, -5], [-100, 1000, 1000, -100, -100], 'k', 'linewidth', 0.1); 
axis off tight; 
% 
if export_fig
    saveas(gcf, sprintf('%s/example_box_trace.fig', output_folder));
    saveas(gcf, sprintf('%s/example_box_trace.pdf', output_folder));
end


% axes('position', [0.1, 0.72, 0.88, 0.27]); 
% plot(t, y1, 'g');
% axis  tight;  
% box on; 
% set(gca, 'ytick', -0:200:400); 
% ylim([-100, 350]); 
% legend('center ROI'); 
% set(gca, 'xticklabel', []); 
% ylabel('Fluo.'); 
% axis off; 
% 
% axes('position', [0.1, 0.42, 0.88, 0.27]); 
% plot(t, y2, 'r'); 
% axis  tight; 
% box on;   
% set(gca, 'ytick', -0:200:400); 
% ylim([-100, 350]); 
% legend('periphery ROI'); 
% set(gca, 'xticklabel', []); 
% ylabel('Fluo.'); 
% axis off; 
% 
% axes('position', [0.1, 0.12, 0.88, 0.27]); 
% plot(t, y1-y2, 'b'); 
% box on; 
% axis  tight; 
% set(gca, 'ytick', -0:200:400); 
% ylim([-100, 350]); 
% legend('difference'); 
% xlabel('Time (sec.)'); 
% ylabel('Fluo.'); 
% axis off; hold on; 
% %% run NMF within the cropped area 
% Ybox = reshape(Ybox, [], T); 
% [u,s,v] = svdsecon(Ybox, 30); 
% 
% figure('papersize', [2, 4]);
% init_fig;
% dw = 1/2; 
% dh = 1/4; 
% 
% for m=1:30
% %     ii = mod(m, 7);
% %     jj = ceil(m/7);
% %     axes('position', [0.5-ii*dw+0.001, 1-jj*dh+0.001, dw-0.001, dh-0.001]);
% 
%     imagesc(reshape(u(:,m), 2*gSiz+1, []));
%     axis equal off tight; 
%     pause; 
% end
% 
% 
% %%
% % temporal 
% figure('papersize', [8, 5]); 
% init_fig; 
% set(gcf, 'defaultAxesFontSize', 16); 
% C = neuron.C_raw(ids_cnmfe_only(1:10),:); 
% C = bsxfun(@times, C, 1./max(C,[],2)); 
% C = bsxfun(@plus, C, (10:-1:1)'); 
% plot((1:T)/neuron.Fs, C','k', 'linewidth', 1); 
% axis tight;
% set(gca, 'ytick', 1:10);
% set(gca, 'yticklabel', 10:-1:1); 
% xlabel('Time (sec.)'); 
% set(gca, 'position', [0.05,0.15,0.92,.84]);
% % ylabel('Neuron #'); 
% box on ; 
% if export_fig
%     saveas(gcf, sprintf('%s/missed_temporal.fig', output_folder) );
%     saveas(gcf, sprintf('%s/missed_temporal.pdf', output_folder));
% end
% 
% %% run rank-1 NMF 
% Ybox = reshape(Ybox, [], T); 
% Ybox_no_baseline = bsxfun(@minus, Ybox, median(Ybox,2)); 
% [A_cnmf, C_cnmf] = nnmf(Ybox_no_baseline, 10); 
% 
% %% display 10 example frames in the cropped box 
% close all; 
% rng(1); 
% Ybox = reshape(Ybox, [], T); 
% Nexample = 5; 
% ind = randi(T, Nexample, 1); 
% w = (2*gSiz+5)*Nexample-5; 
% h = 2*gSiz+1; 
% Ybox_no_baseline = reshape(Ybox_no_baseline, 2*gSiz+1, [], T); 
% figure('papersize', [w, h]/gSiz); 
% init_fig; 
% for m=1:Nexample
%     axes('position', [(2*gSiz+5)*(m-1)/w, 0, (2*gSiz+1)/w, (2*gSiz+1)/h]); 
%     imagesc(Ybox_no_baseline(:, :, ind(m))); 
%     colormap gray; 
%     axis equal off tight; 
% end 
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_box_frames.eps', output_folder), 'psc2');
%     saveas(gcf, sprintf('%s/example_box_frames.pdf', output_folder));
% end
% %% run rank 1 CNMF 
% Ybox_no_baseline = reshape(Ybox_no_baseline, [], T); 
% [W, R] = nnmf(Ybox_no_baseline, 1); 
% % initialized neuron 
% figure('papersize', [4,4]);
% init_fig;
% axes('position', [0, 0, 1, 1]); 
% imagesc(reshape(W, 2*gSiz+1,[])); 
% axis equal tight off;
% hold on;
% 
% if export_fig
%     saveas(gcf, sprintf('%s/example_box_nmf1_spatial.eps', output_folder), 'psc2');
%     saveas(gcf, sprintf('%s/example_box_nmf1_spatial.pdf', output_folder));
% end
% 
% %% show calcium traces 
% figure('papersize', [4.3, 1.2]);
% init_fig;
% axes('position', [0.05, 0.05, 0.9, 0.9]); 
% plot(R, 'b'); 
% axis tight; 
% set(gca, 'xtick', []); 
% set(gca, 'ytick', []); 
% box on; 
% ylim([-0.2, 1.2] *max(R)); 
% if export_fig
%     saveas(gcf, sprintf('%s/example_box_nmf1_temporal.eps', output_folder), 'psc2');
%     saveas(gcf, sprintf('%s/example_box_nmf1_temporal.pdf', output_folder));
% end
% 
% %% 
% 
% %% apply SVD to the raw data
% Y = reshape(Y, [], T); 
% Ymean = mean(Y, 2); 
% [U, S, V] = svdsecon(bsxfun(@minus, Y, Ymean), 30);
% Ymean = mean(Y, 2); 
% % [U, S, V] = svdsecon(Y , 30); 
% 
% %% pick the top 6 singular vectors to show the results 
% 
% Nexample = 5; 
% ind = randi(T, Nexample, 1); 
% w = (2*gSiz+5)*Nexample-5; 
% h = 2*gSiz+1; 
% Ybox_no_baseline = reshape(Ybox_no_baseline, 2*gSiz+1, [], T); 
% figure('papersize', [w, h]/gSiz); 
% init_fig; 
% for m=1:Nexample
%     axes('position', [(2*gSiz+5)*(m-1)/w, 0, (2*gSiz+1)/w, (2*gSiz+1)/h]); 
%     imagesc(Ybox_no_baseline(:, :, ind(m))); 
%     colormap gray; 
%     axis equal off tight; 
% end 
% 
% 
% space_r = 10; 
% space_c = 10; 
% h = d1*2+ space_r*1; 
% w = d2*3 + space_c*3; 
%  
% hf = figure('papersize', [w, h]*0.01); 
% init_fig; 
% for m=1:6
%     [tmp_c, tmp_r] = ind2sub([3,2], m); 
%     axes('parent', hf, 'unit', 'normalized', 'position', [(space_c+d2)*(tmp_c-1)/w, 1+space_r/h-(d1+space_r)*(tmp_r)/h,  d2/w, d1/h]); 
%     neuron.image(sign(sum(U(:,m)))*U(:,m)); 
%     axis equal off tight; 
% end 
% if export_fig
%     saveas(gcf, sprintf('%s/example_top_components_spatial.eps', output_folder), 'psc2');
%     saveas(gcf, sprintf('%s/example_top_components_spatiall.pdf', output_folder));
% end
% 
% %% plot the traces of the top 10 components 
% hf = figure('papersize', [8, 3]); 
% init_fig; 
% axes('position', [0, 0, 1, 1]); 
% hold on; 
% for m=1:6
%     plot(-V(:,m)+m*0.1, 'k', 'linewidth',1 ); 
% end 
% axis tight; 
% % xlim([0, 5000]); 
% axis off; 
% if export_fig
%     saveas(gcf, sprintf('%s/example_top_components_temporal.eps', output_folder), 'psc2'); %#ok<*UNRCH>
%     saveas(gcf, sprintf('%s/example_top_components_temporal.pdf', output_folder));
% end
% 
% %% show SVD 
% figure('papersize', [4,4]); 
% init_fig; 
% s = diag(S); 
% plot(s, '-ob'); 
% xlim([0, 30]); 
% ylim([s(end)*0.8, s(1)*1.5]); 
% set(gca, 'yscale', 'log'); 
% xlabel('Component number'); 
% ylabel('Eigenvalues'); 
% set(gca, 'ytick', [1e4, 1e5, 1e6]); 
% box on; 
% if export_fig
%     saveas(gcf, sprintf('%s/example_eigenvalues.eps', output_folder), 'psc2'); %#ok<*UNRCH>
%     saveas(gcf, sprintf('%s/example_eigenvalues.pdf', output_folder));
% end