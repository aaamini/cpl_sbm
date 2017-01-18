function [priErr, PErr] = compParamErr2(c,e,pri,pih,P,Phat)

K = size(P,1);
n = numel(c);

pri = pri(:);
pih = pih(:);

cm = compCM(c,e,K) / n;

[~,I] = max(cm,[],1);
IDX = sub2ind([K K], I,1:K);
Perm = zeros(K,K);
Perm(IDX) = 1;

priErr = norm(pri - Perm*pih,1) / norm(pri,1);
PErr = norm(P - Perm*Phat*Perm,1) / norm(P,1);
