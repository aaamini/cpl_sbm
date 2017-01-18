function normMuI = compMuI(CM)
% normMUI Computes the normalized mutual information between two clusters
%
% CM   confusion matrix

N = sum(CM(:));
normCM = CM/N; % normalized confusion matrix

IDX = CM == 0; % index of zero elements of CM so that we can avoid them

jointEnt = - sum( (normCM(~IDX)).*log(normCM(~IDX)) );

indpt = sum(normCM,2) * sum(normCM,1);
muInfo = sum(normCM(~IDX) .* log(normCM(~IDX) ./ indpt(~IDX)) );

normMuI = muInfo / jointEnt;