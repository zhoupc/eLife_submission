%% run PCA/ICA
if ~exist('Y', 'var')
    Y = neuron.load_patch_data();
end
tic;
[A_ica, C_ica, time_cost] = run_pca_ica(neuron.reshape(Y, 2), nPCs, nICs, 0.1);
toc;

neuron_ica = Sources2D();
neuron_ica.options = neuron.options;
neuron_ica.A = A_ica;
neuron_ica.C = C_ica;
neuron_ica.C_raw = C_ica;
A_ica_backup = A_ica;
% neuron_ica.viewNeurons();
A_ica_before_trim = neuron_ica.A;
neuron_ica.trimSpatial(0.01, 3);

%% deconv temporal traces 
tic; 
tmp_options = neuron.options.deconv_options; 
tmp_options.optimize_b = false;   % 
neuron_ica.options.deconv_options = tmp_options; 
neuron_ica.deconvTemporal(true);  
fprintf('Time cost in deconvolving temporal traces: %.3f\n', toc); 

%% remove false positives 
neuron_bk = neuron_ica.copy(); 
ind_del = neuron_ica.remove_false_positives(); 
A_ica_before_trim(:, ind_del) = [];
% srt = neuron.orderROIs('snr'); 
% A_ica_before_trim = A_ica_before_trim(:, srt); 
% ind_del = neuron_ica.viewNeurons([], neuron_ica.C_raw);
% A_ica_before_trim(:, ind_del) = [];
