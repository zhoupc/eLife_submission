run_cnmfe_vessel; 

%% show results 
addpath('../functions/');
addpath('~/Dropbox/github/CellSort/'); 
output_folder = '../../Figs/Fig_Vessel_subfigs';
export_fig = false; 


%% run PCA/ICA 
nPCs = 2000; 
nICs = 700;
if ~exist('neuron_ica', 'var')
    tic;
    [A_ica, C_ica] = run_pca_ica(neuron.reshape(Y, 2), nPCs, nICs, 0.1);
    toc;
end
neuron_ica = Sources2D();
neuron_ica.options = neuron.options;
neuron_ica.A = A_ica;
neuron_ica.C = C_ica;
neuron_ica.C_raw = C_ica;
A_ica_backup = A_ica;
neuron_ica.viewNeurons();
neuron_ica_backup = neuron_ica.copy(); 
A_ica_before_trim = neuron_ica.A;
neuron_ica.trimSpatial(0.01, 3);
ind_del = false(size(neuron_ica.A, 2), 1); 
for m=1:size(neuron_ica.A,2) 
    y = neuron_ica.C_raw(m,:);
    [b, sn] = estimate_baseline_noise(y); 
    neuron_ica.C_raw(m,:) = y-b; 
    [c, s] = deconvolveCa(neuron_ica.C_raw(m,:), neuron.options.deconv_options, 'sn', sn); 
    if sum(s>0)==0
        ind_del(m) = true; 
    end 
    neuron_ica.C(m,:) = c; 
end 
neuron_ica.delete(find(ind_del)); 
A_ica_before_trim(:, find(ind_del)) = []; 
neuron_ica.viewNeurons([], neuron_ica.C_raw); 
% %% match neurons between two methods 
% match_cnmfe = pair_neurons(neuron.A, neuron.C, neuron_ica.A, neuron_ica.C); 
% match_ica = pair_neurons(neuron_ica.A, neuron_ica.C, neuron.A, neuron.C); 

%% compute SNRs 
K_cnmfe = size(neuron.C_raw, 1); 
K_ica = size(neuron_ica.C, 1); 

snr_ica = var(neuron_ica.C, 0, 2)./var(neuron_ica.C_raw-neuron_ica.C, 0, 2); 
snr_cnmfe = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2); 

%% sort CNMF-E and PCA/ICA neurons according to its SNR 
[~, srt] = sort(snr_cnmfe, 'descend'); 
neuron.orderROIs(srt); 
snr_cnmfe = snr_cnmfe(srt); 

[~, srt] = sort(snr_ica, 'descend'); 
neuron_ica.orderROIs(srt); 
A_ica_before_trim = A_ica_before_trim(:, srt); 
snr_ica = snr_ica(srt); 

%% match neurons 
match_ica = pair_neurons(neuron_ica.A, neuron_ica.C_raw, neuron.A, neuron.C); 

%% sort neurons 
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

example_neuron = 11;
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

%%
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
Ysignal = neuron.reshape(Ysignal, 1); 

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
neuron.image(Ysignal(:, ind_frame) - neuron.A*neuron.C(:, ind_frame), [-0.5, 0.5]*Yacmax);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
hold on;
plot(coor{1}(1,2:end), coor{1}(2,2:end), 'r', 'linewidth', 2); 
axis equal off tight; colorbar;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_res.fig', output_folder));
    saveas(gcf, sprintf('%s/example_frame_res.pdf', output_folder));
end;

%% variance reduced 
var_0 = var(Y, [], 2);  % variance in the data 
var_nobg = var(Ysignal,[],2); 
var_final = var(Yres, [], 2); 
var_no_signal = var(Y-neuron.A*neuron.C, [],2); 
[u,s,v] = svdsecon(bsxfun(@minus, Ybg, mean(Ybg, 2)), 20); 
var_rank1 = var(Y-u(:,1)*s(1)*v(:,1)', [],2); 
var_rank2 = var(Y-u(:,1:2)*s(1:2,1:2)*v(:,1:2)', [], 2); 
var_rank3 = var(Y-u(:,1:3)*s(1:3,1:3)*v(:,1:3)', [], 2); 
var_rank5 = var(Y-u(:,1:5)*s(1:5,1:5)*v(:,1:5)', [], 2); 
var_rank10 = var(Y-u(:,1:10)*s(1:10,1:10)*v(:,1:10)', [], 2); 

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
breakxaxis([0.2, 0.6]); 
% legend('BG', '1st PC', 'BG without the 1st PC', 'denoised neural signal', 'residual'); 
legend('background', 'denoised neural signal', 'residual'); 

if export_fig
    saveas(gcf, sprintf('%s/variance_explained.fig', output_folder) );
    saveas(gcf, sprintf('%s/variance_explained.pdf', output_folder));
end

%% compute correlation images 
create_correlation_images; 

%%  contours of the detected neurons, CNMF-E
temp = neuron.copy(); 
temp.trimSpatial(0.1, 3); 
Coor_cnmfe = plot_contours(temp.A, Cn_filter, 0.6, 0); 
colormap gray ; 

figure; 
temp = neuron_ica.copy(); 
temp.trimSpatial(0.1, 3); 
Coor_ica  = plot_contours(temp.A,  Cn_filter, 0.6, 0); 
colormap gray; 

%% 
for m=1:length(Coor_ica)
    tmp_contour = Coor_ica{m};
    ind_shift = [2, 0];
    while true
        plot(tmp_contour(1, ind_shift(1):(end-ind_shift(2))), tmp_contour(2,ind_shift(1):(end-ind_shift(2))), '-o');
        title(num2str(m));       
        temp = input('good contours?  (y/n) ', 's');
        if temp=='n'
            ind_shift = input('[left shift, right_shift]   ');
        else
            Coor_ica{m} = tmp_contour(:, ind_shift(1):(end-ind_shift(2)));
            break;
        end
    end
end

%% 
for m=1:length(Coor_ica)
    tmp_contour = Coor_ica{m};
    ind_shift = [2, 0];
    while true
        plot(tmp_contour(1, ind_shift(1):(end-ind_shift(2))), tmp_contour(2,ind_shift(1):(end-ind_shift(2))), '-o');
        title(num2str(m));       
        temp = input('good contours?  (y/n) ', 's');
        if temp=='n'
            ind_shift = input('[left shift, right_shift]   ');
        else
            Coor_ica{m} = tmp_contour(:, ind_shift(1):(end-ind_shift(2)));
            break;
        end
    end
end

%% remove extra pixels
figure;
for m=520:K_cnmfe
%     if ~(ctr(m,1)>129 && ctr(m,2)<128 )
%         continue; 
%     end 
    tmp_contour = Coor_cnmfe{m};
    ind_shift = [2, 0];

    while true
        plot(tmp_contour(1, ind_shift(1):(end-ind_shift(2))), tmp_contour(2,ind_shift(1):(end-ind_shift(2))), '-o');
        title(num2str(m)); 
                temp = input('good contours?  (y/n) ', 's');

        if temp=='n'
            ind_shift = input('[left shift, right_shift]   ');
        else
            Coor_cnmfe{m} = tmp_contour(:, ind_shift(1):(end-ind_shift(2)));
            break;
        end
    end
end

%% missed neurons by PCA/ICA
ind_ica_miss = 1:K_cnmfe; 
ind_ica_miss(match_ica.ind_spatial) = []; 
%% contour plots of CNMF-E and PCA/ICA results. 
figure('papersize', [6, d2/d1*5]); 
init_fig; 
[~, srt] = sort(snr_cnmfe, 'descend'); 
plot_contours(neuron.A(:, srt), Cn_filter, 0.6, 0, [], Coor_cnmfe(srt), 3); 
set(gca, 'position', [0.01, 0.02, 0.9, 0.96]); 
% for m=1:10 
%     temp = Coor_cnmfe{ids_match(m)}; 
%     fill(temp(1,:), temp(2,:), 'b', 'facealpha',1, 'edgecolor', 'none'); 
% end 
% for m=1:10
%     temp = Coor_cnmfe{ind_miss(m)}; 
%     fill(temp(1,:), temp(2,:), 'g', 'facealpha',0.8, 'edgecolor', 'none');
% end 
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
init_fig;
[~, srt] = sort(snr_ica, 'descend'); 
plot_contours(neuron_ica.A(:, srt), Cn_filter, 0.6, 0, [], Coor_ica(srt), 3); 
set(gca, 'position', [0.01, 0.02, 0.9, 0.96]); 
% for m=1:10 
%     temp = Coor_ica{ids_match(m)}; 
%     fill(temp(1,:), temp(2,:), 'b', 'facealpha',0.8, 'edgecolor', 'none'); 
% end 
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
% figure('papersize', [2, 5]);
% init_fig;
% ctr = neuron.estCenter(); 
% % ind_miss = [374,387, 389, 391, 394, 395, 397, 406, 408, 410, 430,...
% %     435, 454, 464, 473,491,476,496,501,510,512,530, 536,541]; 
% ind_miss = [389, 391, 394, 395, 397, 406, 408, 430,...
%     464, 473,491,476,496,536]; 
% 
% for m=1:10
%     axes('position', [0.505-0.5*mod(m,2), 1.003-ceil(m/2)*0.2, 0.49,0.196]);
%     r0 = round(ctr(ind_miss(m), 1));
%     c0 = round(ctr(ind_miss(m), 2));
%     img = neuron.reshape(neuron.A(:, ind_miss(m)),2);
%     if r0<23
%         r0 = 23; 
%     elseif r0>=d1-18
%         r0 = d1-18; 
%     end 
%     if c0<21
%         c0 = 21; 
%     elseif c0>=d2-20 
%         c0 = d2-20; 
%     end 
%     imagesc(img(r0+(-22:18), c0+(-20:20)));
% 
%     axis equal off tight;
%     text(1,6, num2str(m), 'color', 'w', 'fontweight', 'bold', 'fontsize', 15);
% end
figure('papersize', [13,1]);
init_fig;
ctr = neuron.estCenter(); 
% ind_miss = [374,387, 389, 391, 394, 395, 397, 406, 408, 410, 430,...
%     435, 454, 464, 473,491,476,496,501,510,512,530, 536,541]; 
ind_miss = [389, 391, 394, 395, 397, 406, 408, 430,...
    464, 473,491,476,496,536]; 

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
% temp = snr_cnmfe(1:K_match)./snr_ica(1:K_match);  
temp = [9, 10, 14, 16, 17, 21, 23, 27, 28, 34, 38, 41, 59, 63, 66, 69, ...
    70, 72, 79, 83, 89, 119, 124, 104, 109, 102, 101, 165, 223, 97, 165, ...
    167, 169, 174, 210, 217, 251, 262, 276, 280, 282, 296, 307,315,353]; 
temp = sort(temp); 
ids_match = [9, 134, 104, 124, 165, 217, 251, 276, 296, 303, 307]; 
ids_ica = match_cnmfe.ind_max(ids_match); 
ctr = neuron.estCenter(); 
ctr_ica =  neuron_ica.estCenter(); 
%% spatial 
% CNMF-E 
figure('papersize', [2, 0.4]); 
init_fig;
axes('position', [0.025, 0.1, 0.95, 0.85]); 
colorbar('southoutside');  set(gca, 'fontsize', 20); 
axis off; 
saveas(gcf, sprintf('%s/colorbar_cnmfe.fig', output_folder)); 
saveas(gcf, sprintf('%s/colorbar_cnmfe.pdf', output_folder)); 
%% PCA/ICA 
figure('papersize', [2, 0.3]); 
init_fig;
axes('position', [0.05, 0.05, 0.9, 4]);
imagesc([], [-0.2, 1.0]); 
colorbar('southoutside');  set(gca, 'fontsize', 18); 
axis off; 
saveas(gcf, sprintf('%s/colorbar_ica.fig', output_folder)); 
saveas(gcf, sprintf('%s/colorbar_ica.pdf', output_folder)); 

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
% figure('papersize', [18, 6.3]); 
% init_fig; 
% set(gcf, 'defaultAxesFontSize', 16); 
% col = colormap(cool); 
% % CNMF-E
% C = neuron.C_raw(ids_match(1:10),:); 
% % C = bsxfun(@times, C, 1./max(C,[],2)); 
% % C = bsxfun(@plus, C,  (10:(-1):1)'); 
% axes('position', [0.05, 0.15, 0.46, 0.84]); hold on; 
% for m=1:10 
%     y = C(m,:); 
% %     [~, sn] = estimate_baseline_noise(y);  
%     plot((1:T)/neuron.Fs, y/max(y)+10-m,'linewidth', 1, 'color', col(64-4*m,:)); 
%     text(-5, (10-m), num2str(m), 'fontsize', 20); 
% end 
% % ylabel('Neuron #'); 
% % set(gca, 'ytick', 1:10);
% % set(gca, 'yticklabel', 10:-1:1); 
% % xlabel('Time (sec.)'); 
% % box on ; 
% axis off; 
% xlim([0,100]); 
% ylim([-1.2, 10]); 
% % PCA/ICA
% axes('position', [0.52, 0.15, 0.46, 0.84]); 
% C = neuron_ica.C_raw(ids_ica(1:10),:); 
% % C = bsxfun(@times, C, 1./max(C,[],2)); 
% % C = bsxfun(@plus, C, (10:(-1):1)'); 
% hold on ;
% for m=1:10 
%     y = C(m,:); 
% %     [~, sn] = estimate_baseline_noise(y); 
%     plot((1:T)/neuron.Fs, y/max(y)+10-m,'linewidth', 1, 'color', col(64-4*m,:)); 
% end 
% axis tight off;
% % set(gca, 'ytick', 1:10);
% % set(gca, 'yticklabel',[]); 
% % xlabel('Time (sec.)'); 
% % box on ; 
% xlim([0, 100]); 
% ylim([-1.2, 10]); 
% plot([90, 100], [-0.5, -0.5], 'linewidth', 10); 
% text(90, -1.0, '10 sec', 'fontweight', 'bold', 'fontsize', 20); 
% if export_fig
%     saveas(gcf, sprintf('%s/matched_temporal_zoom.fig', output_folder) );
%     saveas(gcf, sprintf('%s/matched_temporal_zoom.pdf', output_folder));
% end

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
%% save video
kt = 4;     % play one frame in every kt frames
save_avi = true;
t_begin = 1; 
t_end = T; 

center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2))*1.25;
range_res = [-1,1]*center_ac;
range_ac = center_ac+range_res;
multi_factor = 8;
center_Y = min(Y(:)) + multi_factor*center_ac;
range_Y = center_Y + range_res*multi_factor;
avi_filename = sprintf('../../Videos/striatum_decompose.avi'); 
cnmfe_save_video; 


%% demixing video 
file_nm = sprintf('../../Videos/blood_vessel_demix.avi'); 
neuron.playAC(file_nm, [], [t_begin, t_end]); 









