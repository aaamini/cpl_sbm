classdef dcBlkMod
% dcBlkMod
    properties
        pri                % prior on community labels
        P                  % edge probability matrix
        n                  % number of nodes
        K                  % number of communities
        lambda             % theoretical expected node degree
        
        thv                % a vector thetas for degree-correction
        thp                % probability of each value in thv
        exTheta            % expected vaue of theta
        
        theta              % the actual theta being used
        c                  % the actual labels being used
        As                 % the actual adjacency matrix being used
    end
    
    methods
        function ob = dcBlkMod(n,K,lambda, LowVal, LowProb)
            ob.K = K;
            ob.lambda = lambda;
            ob.n = n;
            
            ob.pri =  1/K * ones(1,K);
            ob.thv = [1 LowVal];
            ob.thp = [1-LowProb LowProb];
            
            ob.exTheta = ob.thp(:)' * ob.thv(:);
        end
        
        function ob = genP(ob, OIR, inWeights)
            % Generate P matrix with OIR and inWeights
            
            diag_idx = diag(true(ob.K,1));
            
            if OIR == 0
                ob.P = diag(ones(ob.K,1));
            else
                ob.P = ones(ob.K);
                ob.P( diag_idx ) = ob.P(1,2)/OIR;
            end

            ob.P(diag_idx) = inWeights(:).*ob.P(diag_idx);

            %tmplam = (ob.n-1) * ob.pri(:)' * ob.P * ob.pri(:) * (ob.exTheta)^2;
            tmplam = dcBlkMod.compLam(ob.n, ob.pri, ob.P, ob.exTheta);
            ob.P = ob.P *(ob.lambda/tmplam);
            
            if max(ob.P(:)) > 1
                warning('overflow:Pmat','Maximum of P matrix is above 1.')
            end
        end
        
        function ob = genData(ob)
            
            [ob.c,~] = find(mnrnd(1,ob.pri,ob.n)');
            ob.theta = randsrc(ob.n,1,[ob.thv; ob.thp]);
        
            ob.As = genDCBlkMod(ob.c, ob.P, ob.lambda*ob.n, ob.theta);
            ob.As = ob.As + ob.As';
        end
        
        function saveToFile(ob,tag)
            
            if nargin < 2
                tag = '';
            end
            dlmwrite(['graph' tag '.txt'], full(ob.As))
            dlmwrite(['clusters' tag '.txt'], ob.c)
            %dlmwrite(['thetas' tag '.txt'], ob.theta)
          
        end
        
        function ob = perturbData(ob,rho)
            Bs = genBlkMod(ones(ob.n,1),log(ob.n)/ob.n,log(ob.n)*ob.n);
            Bs = Bs + Bs';

            avgDeg = full(mean(sum(ob.As,2))); % average degree of As
            ob.As = ob.As + Bs * rho * avgDeg/log(ob.n);
        end
        
        function thresh = compThresh(ob)
             thresh = zeros(1,ob.K);
             for k = 1:ob.K
                 temp = min(abs( ...
                     repmat(ob.P(k,k), 1, ob.K-1) - ob.P(k,setdiff(1:ob.K,k)) ));
                 thresh(k) = ob.n * temp / (ob.K*sqrt(ob.lambda));
             end
        end
        
        function ob = removeZeroDeg(ob)
            zdNodes = sum(ob.As,2) == 0;
%             nOrg = size(ob.As,1);
            
            ob.As = ob.As(~zdNodes,~zdNodes);
            ob.c = ob.c(~zdNodes); 
            ob.theta = ob.theta(~zdNodes); 
        end
        
        
        function Ph = compPhat(ob)

            K = numel(unique(ob.c));
            Ph = zeros(K);
            IDX = cell(1,K);
            %csizes = zeros(K,1);

            for k = 1:K
                IDX{k} = find(ob.c == k);
            %    csizes(k) = numel(IDX{k});
            end

            for k = 1:K
                for r = 1:K
                    temp  =  ob.As( IDX{k},IDX{r} ) ;
                    Ph(k,r) = mean( temp(:) );
                end
            end
        end
        
        function lamh = compLamh(ob)
            
            nh = size(ob.As,1);
            
            lamh = (1/nh) * sum(ob.As(:));
        end
    end
      
    methods(Static)
        function lam = compLam(n, pri, P, exTheta)
            lam = (n-1) * pri(:)' * P * pri(:) * (exTheta)^2;
        end
    end
end