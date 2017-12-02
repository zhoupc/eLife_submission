Cn = Cn_nobg; 

%% save figures;
T = size(neuron.C,  2);
avi_file = VideoWriter(intervention_video);
avi_file.FrameRate = 1;
avi_file.open();
operation_0 = false(size(neuron_0.A,2),4); 
operation_1 = false(size(neuron_1.A,2),4); 
operation_2 = false(size(neuron_2.A,2),4); 
operation_3 = false(size(neuron_3.A,2),4); 
operation_4 = false(size(neuron_4.A,2),4); 

%% deletion in the first iteration 
ids_before= neuron_0.ids; 
ind = false(size(ids_before)); 
ids_after = neuron_1.ids; 
for m=1:length(ids_before)
    ind(m) = any(ids_before(m)==ids_after); 
end
operation_0(ind==0, 2) = true; 

%% to be merged in the second iteration 
rois_before = [22, 57, 13, 19, 14, 3, 6, 7, 4];
rois_after = [22, 13, 3]; 
ids_before = neuron_1.ids; 
ids_after = neuron_2.ids; 
for m=1:length(rois_before)
    operation_1(ids_before==rois_before(m), 1) = true; 
end
for m=1:length(rois_after)
    operation_2(ids_after==rois_after(m), 3) = true; 
end

%% run more iteration 
neuron_3 = neuron_4; 
% ids_before = neuron_2.ids; 
% ids_after = neuron_4.ids; 
% x = setdiff(ids_before, ids_after); 
% ind = false(size(ids_before)); 
% for m=1:length(x)
%     ind(m) = any(ids_before==x(m));  
% end
% operation_2(ind==0, 1) = true;   % to be deleted 

%% run more iteration 

% %% 
% operation_0([5,7,27], 2) = true; %to be deleted 
% operation_1([6,10,26,27,29,32], 3) = true ; %added 
% operation_1([1,25,29], 4) = true; % merged 
% operation_1([5,6,12,13], 1) = 1; % to be merged  
% operation_1([14,18,26,29,30,32], 2) = 1; %to be deleted
% operation_2([5,12], 4) = true;  %merged 
% operation_2([15,16,19,20], 1) = 1; % to be merged 
% operation_3([15,18], 4) = 1; % merged 
for mm=0:3
    %%
    neuron = eval(sprintf('neuron_%d', mm));
    operation = eval(sprintf('operation_%d', mm));
    K = size(neuron.C, 1);
    % order neurons
    snr = var(neuron.C_raw,0,2)./var(neuron.C_raw-neuron.C, 0, 2);
    [snr, srt] = sort(snr, 'descend');
    neuron.orderROIs(srt);
    % order neurons by letting closeby neurons having similar colors.
    ctr = neuron.estCenter();
    tmp_dist = bsxfun(@minus, ctr(:,1), ctr(:,1)').^2 + bsxfun(@minus, ctr(:,2), ctr(:,2)').^2;
    srt = nan(K,1);
    srt(1) = 1;
    for m=1:(K-1)
        temp = tmp_dist(srt(m), :);
        temp(srt(~isnan(srt))) = inf;
        [~, srt(m+1)] = min(temp);
    end
    neuron.orderROIs(srt);
    
    if exist(sprintf('coor_%d', mm), 'var')
        eval(sprintf('coor = coor_%d;', mm));
    else
        figure;
        coor = plot_contours(neuron.A, Cn, 0.8, 1);
        eval(sprintf('coor_%d=coor;', mm));
        close;
    end
    save_intervention;
    
    axes('position', [0.28, 0.92, 0.69, 0.03]);
    if mm==0
        title('After automated analysis');
%         text(0, 2, 'automated: all', 'fontsize', 16, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
        axes('position', [0.05, 0.63, 0.93, 0.05]); hold on;
        text(25, 6, '{\color{red}to be deleted (manually selected)}', 'fontsize', 15, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
    elseif mm==1
        title('Merge neurons');
%         text(0, 1, 'adding/merging/deleting/post-processing', 'fontsize', 15, 'color', 'k'); 

        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
        axes('position', [0.05, 0.63, 0.93, 0.05]); hold on;
        text(25, 6, '{\color{magenta}to be merged (manually verified)}', 'fontsize', 15, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
    elseif mm==2
    title('After merging neurons');     
%     text(0, 1, 'automated: merging/post-processing', 'fontsize', 15, 'color', 'k'); 
%         text(0, -2, 'manual: deleting', 'fontsize', 15, 'color', 'k'); 
        
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
        axes('position', [0.05, 0.63, 0.93, 0.05]); hold on;
        text(25, 10, '{\color{cyan}merged components}', 'fontsize', 15, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
    elseif mm==3
        title('After one iteration of updating A, C & B');
        %       text(0, 1, 'automated: post-processing', 'fontsize', 15, 'color', 'k');
        %         text(0, -2, 'manual: merging', 'fontsize', 15, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
        axes('position', [0.05, 0.63, 0.93, 0.05]); hold on;
        text(0, 1, '{\color{cyan}merged}', 'fontsize', 15, 'color', 'k');
        ylim([-2, 2]);
        xlim([0, 100]);
        axis off;
    end
    temp = getframe(gcf);
    temp = imresize(temp.cdata, [600, 600]);
    avi_file.writeVideo(temp);
end
avi_file.close();