function C_sim = cos_similarity(C1, C2)
%% compute cosine similarity
C1 = bsxfun(@times, C1, 1./sqrt(sum(C1.^2, 1)));

if nargin<2
    C2 = C1;
else
    C2 = bsxfun(@times, C2, 1./sqrt(sum(C2.^2, 1)));
end

C_sim = C1'*C2;
