%% generate spatial components 
[A, coor] = gen_spatial(K, [d1,d2], gSig, d_min, seed); 
A = neuron.reshape(A, 1); 

%% temporal components 
rng(seed); 

% generate the firing probabilities for each neuron, which is a T*K matrix. P = L*M,
nl = ceil(K/4);     % rank of P 

L = randn(T, nl);    % low dimension states, use random work model 
for m=2:T
    L(m, :) = L(m,:) + L(m-1,:); 
end
L = bsxfun(@minus, L, min(L,[],1)); 
L = bsxfun(@times, L, 1./mean(L,1)/T); 
% generating the mixing matrix 
x0 = (rand(nl,1)-0.5) * 0.8 *d2 + d2/2;    
y0 = (rand(nl,1)-0.5) * 0.8 *d1 + d1/2;
sig = rand(nl,1) * 5 * gSig; 
d_nl = sqrt(bsxfun(@minus, x0, coor.x0').^2 + bsxfun(@minus,y0, coor.y0').^2); 
M = exp(-bsxfun(@times, d_nl, 1./sig).^2); 
M = bsxfun(@times, M, 1./mean(M,1)/nl); 

P = L * M * mu_bar * T; 

% generate 'spike counts' within each frame and simulate its calcium traces
S = zeros(K, T); 
C = zeros(K, T); 
tau_d = tau_d_bar * max(0.5, rand(K,1)*0.2+0.9); 
tau_r = tau_r_bar * max(0.5, rand(K,1)*0.2+0.9); 
for m=1:K 
    % spike counts
    spks = poissrnd(P(:,m)); 
    % merge nearby spikes into 1 to make the spiking signal sparse
    tsp = find(spks);
    t0 = tsp(1);
    for n=2:length(tsp)
        if or(tsp(n)<=t0+minIEI, spks(tsp(n))<minEventSize)
            spks(t0) = spks(t0) + spks(tsp(n));
            spks(tsp(n)) = 0; 
        else
            t0 = tsp(n);
        end
    end
    S(m, :) = spks;

    % simulate spike trains 
    kernel = (exp(-(1:T)/tau_d(m)) - exp(-(1:T)/tau_r(m))) / (tau_d(m)-tau_r(m)); 
    c = conv(spks, kernel); 
    C(m,:) = c(1:T)/std(c(1:T)); 
end 