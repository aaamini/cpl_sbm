function As = genDCBlkMod2(clabels, Pmat, nlambda, theta)
  % generate a degree-corrected block model by probability decomposition

  %disp('Using version 2a')
  
  K = size(Pmat, 1);
  csizes = zeros(1, K);  
  n = length(clabels);
  for i=1:K,
	csizes(i) = sum(clabels == i);
  end

  A1 = KblockGraph(csizes, Pmat);

  % degree-corrected
  [i,j,s]=find(A1);
  ups = (i > j);
  i1 = i(ups);
  j1 = j(ups);
  s1 = s(ups);
  x = s1 .* (rand(length(s1), 1) < theta(i1) .* theta(j1));
  
  
  %[~, Ic] = sort(clabels);
  %Ic_inv = zeros(size(Ic(:)));
  %Ic_inv( Ic ) = 1:n;
  %A2 = sparse([i1;j1],[j1;i1],[x;x],n,n);
  %As = A2(Ic_inv,Ic_inv);

  % adjust symmetry and the order according to clabels
   [~, Ic] = sort(clabels);
   ii = Ic([i1; j1]);
   jj = Ic([j1; i1]);
   As = sparse(ii, jj, [x; x], n, n);
end

function As = KblockGraph(csizes, Pmat)
% K blocks generation
  % csizes: vector of block sizes
  % Pmat: density within/between blocks

  indexrr = [];
  indexcc = [];
  cumcsizes = cumsum(csizes);
  K = length(csizes);
  n = cumcsizes(K);
  for i=1:K,
    for j=i:K,
        
      if (i == j)
        [tmprr, tmpcc] = symblock(csizes(i), Pmat(i, i));
      else
        [tmprr, tmpcc] = offblock(csizes(i), csizes(j), Pmat(i,j));
      end
      
      if (i > 1)
        tmprr = tmprr + cumcsizes(i - 1);
      end
      if (j > 1)
        tmpcc = tmpcc + cumcsizes(j - 1);
      end
      indexrr = [indexrr, tmprr];
      indexcc = [indexcc, tmpcc];
      
    end
  end
  % symmetry
  As = sparse([indexrr, indexcc], [indexcc, indexrr], ones(length(indexrr) * 2, 1), n, n);
end


function [rr, cc] = offblock(nr, nc, p) 
  rr=[]; 
  cc=[];
  if (p > 0) 
      nsize = nr * nc;
      if (nsize < 1e7)
        msample = binornd(nsize, p, 1, 1);
      elseif (nsize*p > 1e6),
        msample = round(normrnd(nsize*p, sqrt(nsize*p*(1-p)),1,1)); 
        msample = max(0, msample);
      else
        msample = poissrnd(nsize*p, 1, 1);
      end
      if (msample > 0)
	    index = mysamplefun(nsize, msample);
        rr = ceil(index/nc);
        cc = index - (rr-1) * nc;
      end
  end
end


function  [rr, cc] = symblock(n, p)
% E-R: (r, c) pairs
rr = [];
cc = [];
if (p>0)
m = n*(n-1)/2;
if (m<1e7),
    msample = binornd(m, p, 1, 1);
elseif (m*p>1e6),
     msample = max(0, round(normrnd(m*p,sqrt(m*p*(1-p)),1,1)));
else
     msample = poissrnd(m*p,1,1);
end 
if (msample > 0)
    index= sort(mysamplefun(m, msample));
    % reorder them
    kk = 1 : (n-1);
    ss = [0.5,  kk .* (kk + 1)./2 + 0.5];
    % icut = as.numeric(cut(index, ss))
    % cc = icut + 1
    cc = ceil(sqrt(1+8*index)/2 + 0.5); 
    rr = index - (ss(cc - 1) - .5);
 end
end 

end


function  vec=mysamplefun(m, n)
% sample n numbers from 1:m
%if (n < 4e4) || (n/m > 0.5)
if (m <= 1e8)
    vec = reshape(randsample(m,n),1,[]);
else
%    error('we should not be here')
    while (1)
    x = ceil(rand(1, n*2) * m);
    y = unique(x);
    if (length(y) >= n)
       break
    end
    end
    %vec = y(1:n);
    vec = randsample(y,n);
end 
end