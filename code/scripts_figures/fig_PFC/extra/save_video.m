y_quantile = 0.9999; 
kt = 10; 
save_avi = true; 
%%
Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
% Y = neuron.reshape(Y, 2);
Ysignal = neuron.reshape(Ysignal, 2);

figure('position', [0,0, 1248, 600], 'color', 'w');
set(gcf, 'defaultAxesFontSize', 12); 
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
if save_avi
    avi_file = VideoWriter('results.avi');
    avi_file.FrameRate = Fs; 
    avi_file.open();
end
temp  = quantile(Ybg(1:1000:numel(Ybg)), [0.0001, y_quantile]);
Ymin = temp(1);
Ymax = temp(2);

ACmax = 300;
ACmin = 0; 

%     subplot(4,6, [5,6,11,12]);
for m=100:kt:T
    subplot(4,6, [1,2, 7, 8]); cla; 
    imagesc(Ybg(:, :,m)+Ysignal(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
        text(d2/2-30, 10, sprintf('Time: %.2f Sec', m/Fs), 'color', 'w');

    subplot(4, 6, [1,2, 7, 8]+12); cla; 
    imagesc(Ybg(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('estimated background');
    colorbar;
    
    subplot(4,6, [3,4, 9,10]); cla; 
    imagesc(Ysignal(:, :, m), [ACmin, ACmax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background-subtracted data'); hold on; colorbar;
    
    
    subplot(4, 6, [3,4, 9,10]+12); cla; 
    imagesc(Yac(:, :, m), [ACmin, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('extracted neural signal'); hold on; colorbar;
    
%         
%     subplot(4,6, [5, 6, 11, 12]);   cla;  
%     neuron.image(neuron_amp* A*C(:, m), [ACmin, ACmax]);
%     axis equal off tight; title('true neural signal'); colorbar;
%     
subplot(4,6, [5,6,11,12]+12); cla; 
    imagesc(Ysignal(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('residual'); hold on; colorbar;
    
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [600, 1248]);
        avi_file.writeVideo(temp);
    end
end

if save_avi
    avi_file.close();
end
