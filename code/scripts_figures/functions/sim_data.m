%% simulate neuron signal 
sim_AC; 

%% generate background sources 
% spatial 
[G, G_coor] = gen_spatial(K_bg, [d1,d2], gSig*bg_neuron_ratio, d_min*5, seed*2); 

% temporal 
F = zeros(K_bg, T); 
for m=1:K_bg
    g = rand(1)*0.3+0.5; 
    tmp_f = randn(1, T); 
    for t=2:T
        tmp_f(t) = tmp_f(t-1)*g + tmp_f(t); 
    end
    tmp_f = tmp_f - min(tmp_f); 
    F(m,:) = tmp_f/std(tmp_f); 
end 

% replace the first background source with a global fluctuations 
load background; 
G(:, :,1) = imresize(G_global, [d1, d2]); 
F(1,:) = resample(F_global, T, length(F_global)); 

% replace the second background source with a blood vessels 
img = zeros(d1,d2); 
x = (1:d2); 
y =  d1 - round(4*(x/d2-1/2).^3*d1 + d1/2); 
y(y<1) = 1; y(y>d1) = d1; 
ind = sub2ind([d1,d2], y, x); 
img(ind) = 1; 
img = imfilter(img, fspecial('gaussian', 10, gSig)); 
G(:, :, 2) = img/max(img(:)); 
F(2,:) = resample(F_vessel, T, length(F_vessel)); 

% change the amplitude of each background source 
ind_center = sub2ind([d1,d2], round(G_coor.y0), round(G_coor.x0)); 
temp = G(:,:,1); 
F = bsxfun(@times, F, 10*[5; 2; 3* sqrt(temp(ind_center(3:end)))]); 
G = neuron.reshape(G,1); 

%% combine results and generate the observed fluorescence 
Fmean = mean(F,2); 
Bf = G * bsxfun(@minus, F, Fmean); 
Bc = G * Fmean * ones(1,T); 
Ysignal = A * C * lambda; 
var_B = sum(Bf.^2, 2)/T; 
rng(seed); 
E = randn(size(Bf)); 
Y = Ysignal + Bf + Bc + sn*E; 












