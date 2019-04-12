function [chat, e] = find_labels(A,K)
    % options for the init method and cpl/upl
    init_opts = struct('verbose',false,'rhoPert',.25);
   
    cpl_opts = struct('verbose',false,'delta_max',0.1, ...
                      'itr_num',20,'em_max',80,'track_err',false);

    %tic, fprintf('%-40s','Applying init. method (SCP) ...')  
    % use spectral clustering with pertubation to initialize labels
    [e init_dT] = initLabel5b(A, K, 'scp', init_opts);    
    
    %tic, fprintf('%-40s','Applying CPL ...') 
    [chat, ~, cpl_dT, post] = cpl4c(A, K, e, [], 'cpl', cpl_opts);
    %chat = e;
end