function [pih, Phat] = cplEstimParam(type, As, post, e)

K = size(post,1);
Phat = zeros(K,K);

switch lower(type)
    case 'hard'
        if any(isnan(post))
            n = size(As,1);
            IDX = sub2ind([K n], e(:)', 1:n);
            post = zeros(K,n);
            post(IDX) = 1;
        end
        pih = mean(post,2);
        pih(isnan(pih)) = 0;
        
        for k = 1:K
            for ell = k:K

                Phat(k,ell) = mean( reshape( As( e == k, e == ell ), [], 1) );
                
                if isnan(Phat(k,ell))
                    Phat(k,ell) = 0;
                end

                 if k ~= ell
                    Phat(ell,k) = Phat(k,ell);
                 end
            end
        end
        
    case 'soft'
          pih = mean(post,2);
          pih(isnan(pih)) = 0;
          
%         n = size(post,2);
%         [~,I]=max(post,[],1);
%         IDX = sub2ind([K n], I,1:n);
%         postNew = zeros(K,n);
%         postNew(IDX) = 1;
%         
        for k = 1:K
            for ell = k:K
                IDXi = post(k,:) >= 0.1;
                IDXj = post(ell,:) >= 0.1;
                
                
                Phat(k,ell) = mean( reshape( ...
                    As(IDXi,IDXj).*(post(k,IDXi)'*post(ell,IDXj)), [], 1) );
                
                if isnan(Phat(k,ell))
                    Phat(k,ell) = 0;
                end

                 if k ~= ell
                    Phat(ell,k) = Phat(k,ell);
                 end
                 
%             disp(sum(sum(post(k,:)'*post(ell,:))))
            end
            
        end
        
        if any(isnan(Phat))
            disp(any(isnan(post)))
            error('y:x','something is wrong')
        end
        
%         postNew
end