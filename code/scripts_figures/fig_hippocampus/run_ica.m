%% run PCA/ICA
if ~exist('Y', 'var')
    Y = neuron.load_patch_data();
end
tic;
rng(10); 
[A_ica, C_ica, time_cost] = run_pca_ica(neuron.reshape(Y, 2), nPCs, nICs, 0.1);
toc;

neuron_ica = Sources2D();
neuron_ica.options = neuron.options;
neuron_ica.A = A_ica;
neuron_ica.C = C_ica;
neuron_ica.C_raw = C_ica;
A_ica_backup = A_ica;
neuron_ica.viewNeurons([]);
neuron_ica_bk = neuron_ica.copy(); 

%% deconv temporal traces 
tic; 
neuron_ica.deconvTemporal(true); 
fprintf('Time cost in deconvolving temporal traces: %.3f\n', toc); 
neuron_ica.orderROIs('srt'); 

%% remove false positives 
A_ica_before_trim = neuron_ica.A;
neuron_ica.options.min_pixel = 25; 
neuron_ica.trimSpatial(0.05, 3); 
neuron_ica.post_process_spatial(); 

ind_del = neuron_ica.remove_false_positives(); 
A_ica_before_trim(:, ind_del) = [];

% neuron_ica.viewNeurons([], neuron_ica.C_raw);

