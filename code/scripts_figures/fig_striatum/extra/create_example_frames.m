figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
Y = neuron.reshape(Y, 1); 
Ybg = neuron.reshape(Ybg, 1); 
Ysignal = neuron.reshape(Ysignal, 1); 

set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
neuron.image(Y(:, ind_frame), [Ymin, Ymax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
colorbar;
axis equal off tight;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame.eps', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_frame.pdf', output_folder));
end

figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ybg(:, ind_frame), [Ymin, Ymax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized'); 
colorbar;
axis equal off tight;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_bg.eps', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_frame_bg.pdf', output_folder));
end

figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(neuron.A*neuron.C(:, ind_frame), [Yacmin, Yacmax]);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized'); 
colorbar;
axis equal off tight;
hold on;

if export_fig
    saveas(gcf, sprintf('%s/example_frame_ac.eps', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_frame_ac.pdf', output_folder));
end

figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
neuron.image(Ysignal(:, ind_frame) - neuron.A*neuron.C(:, ind_frame), [-0.5, 0.5]*Yacmax);
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
axis equal off tight; colorbar;
hold on;
if export_fig
    saveas(gcf, sprintf('%s/example_frame_res.eps', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_frame_res.pdf', output_folder));
end