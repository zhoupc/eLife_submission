%% load previous results 
clear; clc; close all; 
results_prev = matfile('/data/zhoupc/Data/paper_results/pfc_results.mat'); 
Y = results_prev.Y; 
neuron_1 = results_prev.neuron; 
Coor_1 = results_prev.Coor_cnmfe; 
bg_neuron_ratio = 2; 

%% 
C_1 = neuron_1.C; 
A_1 = neuron_1.reshape(neuron_1.A, 2);  
Y_1 = neuron_1.reshape(Y,2); 
[K, T] = size(C_1); 
%% 
for ds_factor = 2:2
    Y = imresize(Y_1, 1.0/ds_factor, 'box'); 
    [tmp_d1, tmp_d2, ~] = size(Y); 
    Y = reshape(Y, [], T);
    neuron = neuron_1.copy(); 
    neuron.P.sn = []; 
    neuron.updateParams('d1', tmp_d1, 'd2', tmp_d2); 
    A = imresize(A_1, 1.0/ds_factor, 'box'); 
    A = reshape(A, [], K);
    neuron.A = A; 
    Y_denoised = A*C_1; 
   
    % estimate background weigths 
    ln = ceil(bg_neuron_ratio*neuron.options.gSiz/ds_factor); 
    sn = imresize(neuron_1.reshape(neuron_1.P.sn, 2), 1/ds_factor, 'box'); 
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
    for m=1:1
        neuron.updateTemporal_endoscope(Ysignal);
    end
    %% save results 
    eval(sprintf('W_%d = W; ', ds_factor)); 
    eval(sprintf('C_%d= neuron.C; ', ds_factor)); 
    eval(sprintf('A_%d=A; ', ds_factor)); 
    
end 

%% compare correlatios 

