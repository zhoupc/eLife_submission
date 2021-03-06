%% correlation image for the raw data
if ~exist('Cn_raw', 'var')
    Y = neuron.reshape(Y, 2);
    Cn_raw = correlation_image(Y, 8);
end
%% correlation image for the background subtracted data
if ~exist('Cn_nobg', 'var')
    Cn_nobg = correlation_image(neuron.reshape(Ysignal, 2), 8);
end
%% correlation image for the spatially filtered data
if ~exist('Cn_filter', 'var')
    [Cn_filter, PNR_filter]  = neuron.correlation_pnr(neuron.reshape(Y, 2));
end

%% plot correlation images
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
imagesc(Cn_raw, [0.85, 1]);
colorbar;
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_Cn_raw.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_Cn_raw.pdf', output_folder));
end

%%
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
imagesc(Cn_nobg, [0, 1]);
colorbar;
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_Cn_nobg.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_Cn_nobg.pdf', output_folder));
end

%%
figure('papersize', [d2+50, d1]/max(d1,d2)*5);
init_fig;
set(gca,'position',[0.01 0.03 .9 0.94],'units','normalized')
imagesc(Cn_filter, [0.6, 1]);
colorbar;
axis equal off tight;

if export_fig
    saveas(gcf, sprintf('%s/example_Cn_filter.fig', output_folder), 'psc2');
    saveas(gcf, sprintf('%s/example_Cn_filter.pdf', output_folder));
end

