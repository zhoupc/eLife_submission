%% data preparation
avi_filename = '../../Videos/bg_comparison_scCNMF_r6.avi'; 
kt = 1; 
t_begin = 1; 
t_end = T; 
save_avi = true; 
center_ac = 40;
range_res = [-1,1]*center_ac;
range_ac = center_ac+range_res;
multi_factor = 5;
center_Y = min(Y(:)) + multi_factor*center_ac;
range_Y = center_Y + range_res*multi_factor;

kernel = strel('disk', 6);
Ysc = neuron.reshape(zeros(size(Y)), 2);
for t=1:T
    Ysc(:, :, t) = imopen(neuron.reshape(Y(:, t), 2), kernel);
end
Ysc = neuron.reshape(Ysc,1); 

%% create avi file
avi_file = VideoWriter(avi_filename);
if ~isnan(neuron.Fs)
    avi_file.FrameRate= neuron.Fs/kt;
end
avi_file.open();

%% play and save
figure('position', [100, 100, 1200, 800]); 
ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
ax_signal_true=   axes('position', [0.015, 0.01, 0.3, 0.42]);
ax_bg=    axes('position', [0.345, 0.51, 0.3, 0.42]);
ax_signal =    axes('position', [0.345, 0.01, 0.3, 0.42]);
ax_bg_nmf =    axes('position', [0.675, 0.51, 0.3, 0.42]);
ax_signal_nmf =     axes('position', [0.675, 0.01, 0.3, 0.42]);
for m=t_begin:kt:t_end
    axes(ax_y); cla; 
    neuron.image(Y(:,m),range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title('Raw data');
    axis equal off tight;
    
    axes(ax_bg); cla; 
    neuron.image(Ybg(:, m),range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
    title('Background, CNMF-E');
    
    axes(ax_signal); cla; 
    neuron.image(Y(:,m)-Ybg(:,m), range_ac); hold on;
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title(sprintf('(Raw-BG) X %d, CNMF-E', multi_factor));
    axis equal off tight;
    
    axes(ax_signal_nmf); cla; 
    neuron.image(Y(:,m)-Ysc(:, m), range_ac);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    title(sprintf('(Raw-BG) X %d, sc-CNMF', multi_factor));
    axis equal off tight;
    
    axes(ax_bg_nmf); cla; 
    neuron.image(Ysc(:, m), range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title('Background, sc-CNMF');
    axis equal off tight;
    %         subplot(4,6, [5,6,11,12]+12);
    
    axes(ax_signal_true); cla;
    neuron.image(Ysignal(:,m), range_ac);  hold on;
    title('Neural signals (ground truth)');
%     text(1, 10, sprintf('Time: %.2f second', m/neuron.Fs), 'color', 'w', 'fontweight', 'bold');
    
    axis equal tight off;
    %     box on; set(gca, 'xtick', []);
    %     set(gca, 'ytick', []);
    drawnow; 
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [800, 1200]);
        avi_file.writeVideo(temp);
    end
end

if save_avi
    avi_file.close();
end
