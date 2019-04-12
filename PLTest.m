function Tstat = PLTest(A, find_labels, K, type, Kp)
    normalized = true;
    
    z = find_labels(A,K);
    if (type == 2)
        y = find_labels(A,Kp);
    else
        Kp=K;
        y=z;
    end
    %y = z;
    Y = label_vec2mat(y);
    b = A * Y;

    Tstat = PLTestStat(b,z,normalized,Kp) ;
end
%%
function Tstat = PLTestStat(b,z,normalized,Kp) 
    d = sum(b,2);
    idx = d > 0;
    n = sum(idx);
    z = z(idx);
    d = d(idx);
    b = b(idx,:);
    Z = label_vec2mat(z);

    Cdeg = Z' * d ;
    if (any( Cdeg == 0)) warning('zero Cdegrees.'), end
    %%
    % remove zero degrees 
    rho = diag(1/Cdeg) * (Z' * b);
    E = diag(d) * (Z * rho);
    idx = E > 0;  
    Tstat = sum( (b(idx) - E(idx)).^2 ./ E(idx)  );
    %%  
    if (normalized) 
        %if (type == 2)
        %    Kb = K;
        %else
            Kb = Kp-1;
        %end
        m = sum(idx);
        C = sqrt(n*Kb);
        Tstat = ((Tstat/C) - C) / sqrt(2);
    end
end 
