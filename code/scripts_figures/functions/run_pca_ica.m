function [A_ica, C_ica, time_cost] = run_pca_ica(Y, nPCs, nICs, mu, clear_workspace)
%% run PCA/ICA method to process calcium imaging data

if ~exist('clear_workspace', 'var') || isempty(clear_workspace)
    clear_workspace = true;
end

%% save data as tiff
fn = 'temp.tif';
if exist(fn, 'file')
    temp = input('use the saved data? [y/n]:   ', 's');
    if ~strcmpi(temp, 'y')
        writeTiff(uint16(Y), fn);
    end
else
    writeTiff(uint16(Y), fn);
end
%% parametes
if ~exist('temp', 'dir')
    mkdir temp;
end
outputdir = './temp';

if ~exist('mu', 'var')
    mu = 0.1; %parameter (between 0 and 1) specifying weight of temporal   information in spatio-temporal ICA
end
termtol = 1e-9;
maxrounds = 1000;

time_cost = zeros(2,1);

%% step 1: perform SVD
tic;
[mixedsig, mixedfilters, CovEvals] = CellsortPCA(fn, [], nPCs, [1,1], outputdir);
time_cost(1) = toc;
fprintf('pca: %.2f\n', time_cost(1));

%% step 2: select which principal components will be kept following dimensional reduction
% PCuse = CellsortChoosePCs(fn, mixedfilters);
PCuse = 1:nPCs;

%% step 3: plot the principal component spectrum and compare with the corresponding random matrix noise floor
% CellsortPlotPCspectrum(fn, CovEvals, PCuse) ;

%% step 4: perform ICA with a standard set of parameters, including skewness as the objective function
tic;
[C_ica, ica_filters] = CellsortICA(mixedsig,mixedfilters, ...
    CovEvals, PCuse, mu, nICs, [], termtol, maxrounds) ;
time_cost(2) = toc;
fprintf('ica: %.2f\n', time_cost(2));
A_ica = reshape(ica_filters, nICs, [])';

if clear_workspace
    delete ./temp/* temp.tif
end



































