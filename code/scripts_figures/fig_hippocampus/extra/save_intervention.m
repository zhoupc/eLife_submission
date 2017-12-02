[K, t] = size(neuron.C); 
tt = (1:T)/neuron.Fs; 
%% 
figure('papersize', [10, 8]); 
init_fig; 
set(gcf, 'defaultaxesfontsize', 16); 
axes('position', [0.05, 0.75, 0.2, 0.2]); colormap(gca, winter); 
coor = neuron.get_contours(0.8); 
plot_contours(neuron.A, Cn, 0.8, 1, [], coor, 2); 
axis equal off;
axis([0.5, d2+0.5, 0.5, d1+0.5]); 
title('Contours'); 
% axes('position', [0.23, 0.65, 0.2, 0.3]); 
% imagesc(Cn_res, [0, 1]); 
% axis equal off tight; 
% title('Corr. image'); 

% neuron shapes 
ctr = neuron.estCenter();
n_per_row = 13;
x = 0:31; 
y = x; 
dw = 0.70/n_per_row; 
dh = 0.2/3; 
l = round(gSiz*1.5); 
for m=(1:K)
    x0 = round(ctr(m, 2)); 
    y0 = round(ctr(m,1)); 
    c0 = min(max(1, x0-l), d2-2*l-1); 
    r0 = min(max(1, y0-l), d1-2*l-1); 
    img = neuron.reshape(neuron.A(:, m), 2); 
    img = img(r0+(1:(2*l+1)), c0+(1:(2*l+1)));  
    jj = mod(m, n_per_row); 
    if jj==0
        jj = n_per_row; 
    end 
%     axes('position', [0.02, 0.05, 0.45, 0.45]); 

    ii = (m-jj)/n_per_row;
    axes('position', [dw*(jj-1)+0.28, dh*(3-ii)+0.68, dw-0.001, dh-0.001]); 
    
    imagesc(img); 
    text(1, 5, num2str(m), 'color', 'w', 'fontsize', 10); hold on; 
    if operation(m, 1) % to be merged 
        plot([1,2*l+1, 2*l+1, 1, 1], [1,1,2*l+1, 2*l+1, 1], 'm'); 
    elseif operation(m,2) % to be deletd 
        plot([1,2*l+1, 2*l+1, 1, 1], [1,1,2*l+1, 2*l+1, 1], 'r'); 
    elseif operation(m,3) % added 
        plot([2,2*l, 2*l, 2, 2], [2,2,2*l, 2*l, 2], 'c'); 
    elseif operation(m,4) % merged 
        plot([2,2*l, 2*l, 2, 2], [2,2,2*l, 2*l, 2], 'g'); 
    end
    axis equal off; 
    axis([0, 2*l+2, 0, 2*l+2]); 
end 
%traces 
axes('position', [0.05, 0.08, 0.93, 0.62]); hold on;
temp = bsxfun(@times, neuron.C_raw, 1.5./max(neuron.C_raw, [], 2)); 
K = size(temp,1); 
temp = bsxfun(@plus, temp, (1:K)'); 
cmap = jet(K);
for m=1:K
    plot(tt, temp(m, :), 'color', cmap(m, :)); 
end 
set(gca, 'ytick', 1:2:K); 
box on; 
set(gca, 'fontsize', 16); 
xlabel('Time (sec.)'); 
axis([0, 610, 0.5, K+1.5]); 