function M = compCM(c,e,K)
% Compute the confusion matrix between labels "c" and "e"
%
% c,e Two sets of labels
% K   number of labels in both "c" and "e"

M = zeros(K);
for k = 1:K
    for r = 1:K
        M(k,r) = sum( (c(:) == k) .* ( e(:) == r ) );
    end
end