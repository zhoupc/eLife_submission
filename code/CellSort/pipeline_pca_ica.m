%% step 0: choose data 
fn = 'sim_data.tif'; 
flims = []; 
nPCs = 275; 
dsamp = [1, 1]; 
outputdir = '~/Dropbox/github/CellSort/'; 
mu = 0.1; %parameter (between 0 and 1) specifying weight of temporal   information in spatio-temporal ICA
nIC = 250; 
termtol = 1e-9; 
maxrounds = 200; 
mode = 'series'; 
tlims = [1, 2000]; 
dt = 1; 
ratebin = 1; 
plottype = 1; 
ICuse = 1:nIC; 
smwidth = 4; % standard deviation of Gaussian smoothing kernel 
thresh = 1; 
arealims = [15, 15]; 
plotting = 0; 
subtractmean = true; 

%% step 1: perform SVD 
[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn,...
 flims, nPCs, dsamp, outputdir);
 
%% step 2: select which principal components will be kept following dimensional reduction 
% PCuse = CellsortChoosePCs(fn, mixedfilters); 
PCuse = 5:nPCs; 
%% step 3: plot the principal component spectrum and compare with the corresponding random matrix noise floor 
CellsortPlotPCspectrum(fn, CovEvals, PCuse) ; 

%% step 4: perform ICA with a standard set of parameters, including skewness as the objective function 
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig,mixedfilters, ...
CovEvals, PCuse, mu, nIC, [], termtol, maxrounds) ; 

%% step 5: cell_ica plot 
ICuse = 11:20; 
CellsortICAplot(mode, ica_filters, ica_sig, movm, tlims, dt, ratebin, plottype, ICuse); 
pause; 

%% step 6: segment spatial filters derived by ICA 
 [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, ...
     smwidth, thresh, arealims, plotting);  
 
 %% step 7: read in movie data and output signals corresponding to specified spatial 
%  cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean); 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 