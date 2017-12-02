
%% create video showing overlapping neurons 
Y = neuron.reshape(Y,2); 
Ysignal = neuron.reshape(Ysignal,2); 
xmin = round(xmin); 
xmax = round(xmax); 
ymin = round(ymin); 
ymax = round(ymax); 
nr = ymax-ymin+1; 
nc = xmax-xmin+1; 
Ybox = Y(ymin:ymax, xmin:xmax,:); 
Ysignal_box = Ysignal(ymin:ymax, xmin:xmax,:); 
temp = neuron.reshape(neuron.A*neuron.C,2); 

Coor_cnmfe_box = Coor_cnmfe(ids); 
for m=1:length(ids)
    Coor_cnmfe_box{m} = bsxfun(@minus, Coor_cnmfe_box{m}, [xmin; ymin]); 
end 
Coor_ica_box = Coor_ica(ica_ids); 
for m=1:length(ica_ids)
    Coor_ica_box{m} = bsxfun(@minus, Coor_ica_box{m}, [xmin; ymin]); 
end 
% spatial components within a box 
temp = neuron.reshape(neuron.A(:, ids),2); 
Abox_cnmfe = reshape(temp(ymin:ymax, xmin:xmax, :), nr*nc, []); 
temp = neuron.reshape(A_ica_before_trim(:, ica_ids),2); 
Abox_ica = reshape(temp(ymin:ymax, xmin:xmax, :), nr*nc, []); 

% colormap 
col_cnmfe = [1,0,0;0,1,0;0,0,1]; 
col_ica = [1,0,0;0,1,0;0,0,1]; 

Y_mixed_cnmfe = zeros(nr*nc, T, 3);
Y_mixed_ica = zeros(nr*nc, T, 3);

for m=1:3
    Y_mixed_cnmfe(:, :, m) = Abox_cnmfe* (diag(col_cnmfe(:,m))*neuron.C(ids,:));
    Y_mixed_ica(:, :, m) = Abox_ica* (diag(col_ica(:,m))*neuron_ica.C_raw(ica_ids,:));
end

center_ac = min(max(neuron.A(:, ids),[],1)'.*max(neuron.C(ids, :),[],2))/3;
range_res = [-1,1]*center_ac;
range_ac = center_ac+range_res;
multi_factor = 6;
center_Y = min(Ybox(:)) + multi_factor*center_ac;
range_Y = center_Y + range_res*multi_factor;

Y_mixed_cnmfe = int8(Y_mixed_cnmfe); 
Y_mixed_ica = int16(Y_mixed_ica); 
%% 
kt = 5; 
t_begin = 1; 
t_end = 9000; 
save_avi = true; 
avi_filename = '../../Videos/PFC_overlapping.avi';
avi_file = VideoWriter(avi_filename);
if ~isnan(neuron.Fs)
    avi_file.FrameRate= neuron.Fs/kt;
end
avi_file.open();

figure('position', [500, 500, 800, 400]); colormap gray; 
ax_raw = axes('position', [0.01, 0.55, 0.22, 0.35]); 
ax_signal = axes('position', [0.01, 0.1, 0.22, 0.35]); 
ax_cnmfe = axes('position', [0.26, 0.55, 0.22, 0.35]); 
ax_ica = axes('position', [0.26, 0.1, 0.22, 0.35]); 
ax_cnmfe_trace = axes('position', [0.51, 0.55, 0.5, 0.35]); 
ax_ica_trace = axes('position', [0.511, 0.1, 0.5, 0.35]); 

% show raw data 
axes(ax_raw); 
temp_y = imagesc(Ybox(:, :, 1), range_Y); hold on; 
plot(Coor_cnmfe_box{1}(1,:), Coor_cnmfe_box{1}(2,:), 'r', 'linewidth', 2); 
plot(Coor_cnmfe_box{2}(1,:), Coor_cnmfe_box{2}(2,:), 'g', 'linewidth', 2); 
plot(Coor_cnmfe_box{3}(1,:), Coor_cnmfe_box{3}(2,:), 'b', 'linewidth', 2); 
axis equal off tight; 
title('Raw data'); 
% show bg-subtraced 
axes(ax_signal); 
temp_signal = imagesc(Ysignal_box(:, :, 1), range_ac); hold on; 
plot(Coor_cnmfe_box{1}(1,:), Coor_cnmfe_box{1}(2,:), 'r', 'linewidth', 2); 
plot(Coor_cnmfe_box{2}(1,:), Coor_cnmfe_box{2}(2,:), 'g', 'linewidth', 2); 
plot(Coor_cnmfe_box{3}(1,:), Coor_cnmfe_box{3}(2,:), 'b', 'linewidth', 2); 
axis equal off tight;
title('(Raw-BG) X 4'); 

% show CNMF-E data 
axes(ax_cnmfe); 
x = reshape(Y_mixed_cnmfe(:, 1, :), nr, nc, 3); 
temp_cnmfe = imagesc(x); hold on; 
plot(Coor_cnmfe_box{1}(1,:), Coor_cnmfe_box{1}(2,:), 'r', 'linewidth', 2); 
plot(Coor_cnmfe_box{2}(1,:), Coor_cnmfe_box{2}(2,:), 'g', 'linewidth', 2); 
plot(Coor_cnmfe_box{3}(1,:), Coor_cnmfe_box{3}(2,:), 'b', 'linewidth', 2); 
axis equal tight; 
title('CNMF-E');
set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
box on; 
% show ICA data 
axes(ax_ica); 
x = reshape(Y_mixed_ica(:, 1, :), nr, nc, 3); 
temp_ica = imagesc(x); hold on; 
plot(Coor_cnmfe_box{1}(1,:), Coor_cnmfe_box{1}(2,:), 'r', 'linewidth', 2); 
plot(Coor_cnmfe_box{2}(1,:), Coor_cnmfe_box{2}(2,:), 'g', 'linewidth', 2); 
plot(Coor_cnmfe_box{3}(1,:), Coor_cnmfe_box{3}(2,:), 'b', 'linewidth', 2); 
axis equal  tight; 
title('PCA/ICA'); 
set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
box on; 
% show CNMF-E traces 
axes(ax_ica_trace); hold on; 
y1 = neuron_ica.C_raw(ica_ids(1), :);
plot((1:T)/neuron.Fs,y1/max(y1)+1, 'r', 'linewidth', 2);
y1 = neuron_ica.C_raw(ica_ids(2), :);
plot((1:T)/neuron.Fs,y1/max(y1)+2, 'g', 'linewidth', 2);
y1 = neuron_ica.C_raw(ica_ids(3), :);
plot((1:T)/neuron.Fs,y1/max(y1)+3, 'b', 'linewidth', 2);
t_cnmfe = plot([t_begin, t_begin], [0.5, 4], 'k');
xlim([-100, 100]);
ylim([0.5, 4]); 
xlabel('Time (sec.)'); 
% show CNMF-E traces
axes(ax_cnmfe_trace); hold on;
y1 = neuron.C_raw(ids(1), :);
plot((1:T)/neuron.Fs,y1/max(y1)+1, 'r', 'linewidth', 2);
y1 = neuron.C_raw(ids(2), :);
plot((1:T)/neuron.Fs,y1/max(y1)+2, 'g', 'linewidth', 2);
y1 = neuron.C_raw(ids(3), :);
plot((1:T)/neuron.Fs,y1/max(y1)+3, 'b', 'linewidth', 2);
t_ica = plot([t_begin, t_begin], [0.5, 4], 'k');
xlim([-100, 100]);
ylim([0.5, 4]); 
for m=t_begin:kt:t_end
    % show raw trace
    axes(ax_raw); delete(temp_y);
    temp_y = imagesc(Ybox(:, :, m), range_Y);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
    
    % show bg-subtraced
    axes(ax_signal); delete(temp_signal);
    temp_signal = imagesc(Ysignal_box(:, :, m), range_ac);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
    pause(0.1);
    
    % show CNMF-E video
    axes(ax_cnmfe); delete(temp_cnmfe);
    x = reshape(Y_mixed_cnmfe(:, m, :), nr, nc, 3);
    temp_cnmfe = imagesc(x);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
   
    % show PCA/ICA video
    axes(ax_ica); delete(temp_ica);
    x = reshape(Y_mixed_ica(:, m, :), nr, nc, 3);
    temp_ica = imagesc(x);     
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;

    % shift time
    axes(ax_cnmfe_trace); 
    delete(t_cnmfe);
    t_cnmfe = plot([m, m]/neuron.Fs, [0.5, 4], 'k');
    xlim([-100, 100]+m/neuron.Fs);
    set(gca, 'ytick', [1,2,3]); 
    
    axes(ax_ica_trace); 
    delete(t_ica);
    t_ica = plot([m, m]/neuron.Fs, [0.5, 4], 'k');
    xlim([-100, 100]+m/neuron.Fs);
    pause(0.1);
        set(gca, 'ytick', [1,2,3]); 

    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [400, 800]);
        avi_file.writeVideo(temp);
    end
end 
avi_file.close(); 