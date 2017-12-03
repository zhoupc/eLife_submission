%% load previous results 
clear; clc; close all; 
results_prev = matfile('/data/zhoupc/Data/paper_results/vessel_results.mat'); 
Y = results_prev.Y; 
neuron_0 = results_prev.neuron; 
Coor_0 = results_prev.Coor_cnmfe; 
bg_neuron_ratio = 2; 

%% 
C_0 = neuron_0.C; 
A_0 = neuron_0.reshape(neuron_0.A, 2);  
Y_0 = neuron_0.reshape(Y,2); 
[K, T] = size(C_0); 
%% 
for ds_factor = 1:8
    Y = imresize(Y_0, 1.0/ds_factor, 'box'); 
    [tmp_d1, tmp_d2, ~] = size(Y); 
    Y = reshape(Y, [], T);
    neuron = neuron_0.copy(); 
    neuron.P.sn = []; 
    neuron.updateParams('d1', tmp_d1, 'd2', tmp_d2); 
    neuron.options.deconv_options = 'constrained'; 
    A = imresize(A_0, 1.0/ds_factor, 'box'); 
    A = reshape(A, [], K);
    neuron.A = A; 
    Y_denoised = A*C_0; 
   
    % estimate background weigths 
    ln = ceil(bg_neuron_ratio*neuron.options.gSiz/ds_factor); 
    sn = imresize(neuron_0.reshape(neuron_0.P.sn, 2), 1/ds_factor, 'box'); 
    Y_bg = reshape(Y-Y_denoised, [tmp_d1, tmp_d2, T]); 
    [Yest, results] = local_background(Y_bg, 1, ln, [], sn, 10); 
    b0 = mean(neuron.reshape(Yest, 1),2); 
    % create a sparse matrix W given the results of background estimation 
    tmp_d = tmp_d1*tmp_d2; 
    nmax = size(results.weights{floor(tmp_d1/2), floor(tmp_d2/2)}, 2); 
    ii = zeros(tmp_d*nmax, 1); 
    jj = ii; 
    ss = ii; 
    k = 0; 
    for m=1:tmp_d
        temp = results.weights{m}; 
        tmp_k = size(temp, 2); 
        ii(k+(1:tmp_k)) = m; 
        jj(k+(1:tmp_k)) = temp(1,:); 
        ss(k+(1:tmp_k)) = temp(2,:); 
        k = k+tmp_k; 
    end 
    ii((k+1):end) = []; 
    jj((k+1):end) = []; 
    ss((k+1):end) = []; 
    W = sparse(ii, jj, ss, tmp_d, tmp_d); 
    %% estimate C given A and W 
    % initialize results 
    Yp = Y - W*Y; 
    Ap = A - W*A; 
    bp = b0 - W*b0; 
    Ysignal = (Yp-bp*ones(1,T)); 
    neuron.C = (Ap'*Ap)\(Ap'*Ysignal); 
    neuron.C(neuron.C<0) = 0; 
    neuron.A = A; 
    Ysignal = Yp - bp*ones(1,T) - (W*A)*neuron.C; 
    % udpate temporal traces 
    neuron.A = Ap; 
    for m=1:1
        neuron.updateTemporal_endoscope(Ysignal);
    end
    %% save results 
    neuron.A = A; 
    neuron.C = bsxfun(@times, neuron.C, neuron.P.sn_neuron'); 
    neuron.C_raw = bsxfun(@times, neuron.C_raw, neuron.P.sn_neuron'); 
    neuron.S = bsxfun(@times, neuron.S, neuron.P.sn_neuron'); 
    eval(sprintf('W_%d = W; ', ds_factor)); 
    eval(sprintf('neuron_%d= neuron.copy(); ', ds_factor)); 
    eval(sprintf('A_%d=A; ', ds_factor)); 
%     
%     %% run 1 phase imaging 
%         fprintf('decimation factor: %d\n', ds_factor);
% 
%     neuron_1phase = neuron.copy(); 
%     neuron_1phase.options.min_pixel = 1; 
%     % initialization 
%     for m=1:K
%         ai = neuron.reshape(A(:, m), 2); 
%         ai_dilate = imfilter(ai, ones(3));
%         ind_cell = (ai>0); 
%         ind_bg = (ai_dilate>0) & (~ind_cell);
%         y_bg = mean(Y(ind_bg(:), :)); 
%         y_cell = mean(Y(ind_cell(:), :))-y_bg; 
%         % deconv y_cell
%         [b, tmp_sn] = estimate_baseline_noise(y_cell);
%         y_cell = y_cell -b;
%         [ck, sk, deconv_options]= deconvolveCa(y_cell, neuron_1phase.options.deconv_options, 'sn', tmp_sn);
%         neuron_1phase.C_raw(m, :) = y_cell; 
%         neuron_1phase.C(m,:) = ck; 
%         
%         Y_box = Y(ind_cell(:), :);
%         X = [y_cell', y_bg', ones(T,1)];
%         temp = (X'*X)\(X'*Y_box');
%         neuron_1phase.A(ind_cell(:),m) = temp(1,:);
%     end
%     
%     [Ybg, results] = local_background(neuron.reshape(Y-neuron_1phase.A*neuron_1phase.C, 2), 1, ln, [], [], 10);
%     Ysignal = Y - neuron.reshape(Ybg, 1);
%     for m=1:2
%         %spatial
%         neuron_1phase.updateSpatial_endoscope(Ysignal, 10, 'hals');
%         
%         %temporal
%         neuron_1phase.updateTemporal_endoscope(Ysignal);
%     end
%     eval(sprintf('neuron_%d_1phase= neuron_1phase.copy(); ', ds_factor));
%     
    %%
end

save ds_results neuron_* 
%% compare correlatios 
corr_denoised = zeros(K,8); 
corr_raw = zeros(K,8); 
corr_deconvolved = zeros(K,8);
for m=1:8
    eval(sprintf('neuron = neuron_%d;', m));
    corr_denoised(:, m) = (mean(neuron_1.C.*neuron.C,2)-mean(neuron_1.C,2).*mean(neuron.C,2))./std(neuron_1.C,0,2)./std(neuron.C,0,2);
    corr_raw(:, m) = (mean(neuron_1.C_raw.*neuron.C_raw,2)-mean(neuron_1.C_raw,2).*mean(neuron.C_raw,2))./std(neuron_1.C_raw,0,2)./std(neuron.C_raw,0,2);
    corr_deconvolved(:, m) = (mean(neuron_1.S.*neuron.S,2)-mean(neuron_1.S,2).*mean(neuron.S,2))./std(neuron_1.S,0,2)./std(neuron.S,0,2);
end

%% 
col = cool(K);
figure('papersize', [8, 8]); 
init_fig; hold on; 
for m=1:K 
    plot(1:8, corr_denoised(m, :), 'color', col(K-m+1, :)); 
end 
errorbar(1:8, mean(corr_denoised, 1), std(corr_denoised, 0, 1), 'r', 'linewidth', 3); 
axis tight; box on; 
set(gca, 'xtick' , 1:8); 
set(gca, 'xticklabel', {'1x1', '2x2', '3x3', '4x4', '5x5', '6x6', '7x7', '8x8'}); 
saveas(gcf, 'corr_denoised_ds_factor.pdf'); 
%% 
id = 1; 
figure('papersize', [8.5, 1.5]);
init_fig; 
for ds_factor=1:8
    axes('position', [(ds_factor-1)*0.12+0.001, 0.001, 0.12, 0.8]); 
    eval(sprintf('neuron = neuron_%d;', ds_factor));
    img = neuron.reshape(neuron.A(:,id),2); 
    neuron.image(img);  axis  off equal ;
    [r0, c0, ~] = find(img==max(img(:))); 
    xlim([c0+[-1,1]*neuron.options.gSiz/ds_factor]); 
    ylim([r0+[-1,1]*neuron.options.gSiz/ds_factor]); 
    title(sprintf('%dx%d', ds_factor, ds_factor)); 
   
end 
saveas(gcf, 'example_neuron.pdf'); 

%% example frame 
 1; 
figure('papersize', [8, 1.41]);
init_fig; 
for ds_factor=1:8
    axes('position', [(ds_factor-1)*0.125+0.001, 0.001, 0.123, 0.7]); 
    eval(sprintf('neuron = neuron_%d;', ds_factor));
    img = Y_0(:, :, 100);
    img = imresize(img, 1.0/ds_factor); 
    neuron.image(img);  axis  off equal tight;
    title(sprintf('%dx%d', ds_factor, ds_factor)); 
   
end 
saveas(gcf, 'example_frame.pdf'); 
%% 
ids = 1;
ds_factors = [2,4, 8]; 
t = (1:T)/neuron.Fs; 
for m_ids=ids
    figure('papersize', [8, 3]); 
    init_fig; 
    axes('position', [0.1, 0.1, 0.9, 0.85]); 
    plot(t, neuron_1.C(m_ids, :)); 
    hold on; 
    cmax = max(neuron_1.C(m_ids,:)); 
    text(10, cmax*1.05,  '1x1, 1.0000'); 
    for m = 1:length(ds_factors)
        m_ds = ds_factors(m); 
        eval(sprintf('neuron = neuron_%d.copy();', m_ds)); 
        plot(t, neuron.C(m_ids, :), '-.'); 
        text(10+m*T/neuron.Fs/5, cmax*1.05, sprintf('%dx%d, %.4f', m_ds, m_ds, corr_denoised(m_ids, m_ds))); 
    end 
%     legend('show', 'orientation', 'horizental'); 
end 























