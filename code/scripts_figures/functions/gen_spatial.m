function [A, coor] = gen_spatial(K, siz, gSig, d_min, seed)
%% simulate spatial components of K neurons within the field of view (siz(1) x siz(2) pixels).
%% input:
%   siz: [d1,d2,T] vector, dimension of the data
%   K: scalar, number of neurons
%   gSig: sigma for the 2D gaussian shape
%   d_min: scalar, minimum distance between two neurons
%   seed: scalar, seed number for generating random numbers

%% output:
%   A: d1*d2*K 3D matrix, spatial component
%   coor: K*2 matrix, center location of each neuron

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

%% initialize parameters
if ~exist('d_min', 'var') || isempty(d_min)
    d_min = 5;
end
if ~exist('seed', 'var')||isempty(seed)
    seed = 1;
end
d1 = siz(1);
d2 = siz(2);
rng(seed);

%% generate spatial profiles for each neuron
% locations of the centers.
x0 = (rand(3*K,1)-0.5) * 0.8 *d2 + d2/2;  % some positions should be removed. thus initialize only 3K
y0 = (rand(3*K,1)-0.5) * 0.8 *d1 + d1/2;

while true
    % generate distances between neurons
    d_neurons = sqrt(bsxfun(@minus, x0, x0').^2+bsxfun(@minus, y0, y0').^2);
    [indr, indc] = find(d_neurons<d_min);
    ind = unique(indr(indr<indc));
    
    % remove pixels that are too close
    if isempty(ind)
        if length(x0)<K
            x0 = [x0; (rand(3*K,1)-0.5) * 0.8 *d2 + d2/2];  %#ok<*AGROW> % some positions should be removed. thus initialize only 3K
            y0 = [y0; (rand(3*K,1)-0.5) * 0.8 *d1 + d1/2];
        else
            break;
        end
    else
        x0(ind) = [];
        y0(ind) = [];
    end
end

x0 = x0(1:K);   % select only top K neurons
y0 = y0(1:K);

% choose sig_x and sig_y for each neuron
sig_x = (1+randn(K,1)*0.1) * gSig;
sig_y = (1+randn(K,1)*0.1) * gSig;
sig_x(or(sig_x<0.5*gSig, sig_x>2*gSig)) = 1;
sig_y(or(sig_y<0.5*gSig, sig_y>2*gSig)) = 1;
[x, y] = meshgrid(1:d2, 1:d1);
A = zeros(d1,d2,K);
coor.x0 = round(x0);
coor.y0 = round(y0);
for m=1:K
    % initialize each neuron's spatial component
    tmp1 = ((x-x0(m))/sig_x(m)).^2 + ((y-y0(m))/sig_y(m)).^2;
    temp = exp(-tmp1/2);
    temp(temp<exp(-5)) = 0;    % remove boundaries
    A(:, :, m) = temp/max(temp(:));
end
