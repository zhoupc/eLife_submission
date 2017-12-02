Y = neuron.reshape(Y, 2);
Ysiz = size(Y);
nam = sprintf('%sdata_snr_1.mat', results_folder);
if ~exist('nam', 'file')
    save(nam, 'Y', 'Ysiz', '-v7.3');
end
neuron.select_data(nam);
Y = neuron.reshape(Y,1);
%% run CNMF-E
rng(1);
neuron.getReady(pars_envs);
[center, Cn, PNR] = neuron.initComponents_parallel(K, [], save_initialization, use_parallel, use_prev);
neuron.merge_neurons_dist_corr();

%% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();
RSS_init = neuron.compute_RSS();
disp(RSS_init); 

%% pick neurons from the residual
neuron.initComponents_residual_parallel([], save_initialization, use_parallel);
neuron.merge_neurons_dist_corr();
RSS_residual = neuron.compute_RSS();
disp(RSS_residual); 

%% correct the inaccurate extraction of neural traces
use_c_hat = false;
neuron.update_temporal_parallel(use_parallel, use_c_hat);
neuron.update_spatial_parallel(use_parallel);
neuron.merge_neurons_dist_corr();
RSS_iter0 = neuron.compute_RSS();

%% udpate spatial&temporal components,
for m=1:8
    disp(m);
    if mod(m, 3) ==0
        neuron.update_background_parallel(use_parallel);
    else
        % update temporal
        neuron.update_temporal_parallel(use_parallel);
        % update _spatial
        neuron.update_spatial_parallel(use_parallel);
    end
    tmp_RSS = neuron.compute_RSS();
    eval(sprintf('RSS_iter_%d=tmp_RSS;', m));
    
    disp(tmp_RSS);
end

RSS_cnmfe = [RSS_init, RSS_residual, RSS_iter0, RSS_iter_1, RSS_iter_2, RSS_iter_3, ...
    RSS_iter_4, RSS_iter_5, RSS_iter_6, RSS_iter_7, RSS_iter_8];
RSS_ids = 1:length(RSS_cnmfe);
labels = [1, 2, 3, 4, 4, 5, 4, 4, 5, 4, 4];
label_meanings = {'initialize A, C, W & b_0', 'add neurons from the residual', ...
 'fix the correlation issue', 'update A, C & b_0 ', 'update W & b_0 '};

%%
% save CNMFE_RESULTS neuron_init neuron D F  A C sn -append;
neuron0_init = neuron_init.copy();
neuron0 = neuron.copy();
neuron0_init.compress_results();
neuron0.compress_results();
save(results_file, 'neuron0_init', 'neuron0', 'D', 'F', 'A', 'C', 'sn', 'RSS_cnmfe', 'RSS_ids', 'labels', 'label_meanings', 'RSS_ground_truth', 'var_Y',  '-append');
