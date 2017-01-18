tic, fprintf('%-40s','Setting up the model ...')
n = 8600; % number of nodes
K = 3;    % number of communities
oir = 0.05;    % "O"ut-"I"n-"R"atio 
lambda = 7;    % average node degree

compErr = @(c,e) compMuI(compCM(c,e,K));    % use mutual info as a measure of error/sim.

inWei = [1 5 10];   % relative wieght of in-class probabilities
[lowVal lowProb] = deal(0.2,0.9); % An example of Degree corrected block model
%[lowVal lowProb] = deal(0.2,0); % An example of original block model


mo = dcBlkMod2(n,K,lambda, lowVal, lowProb); % create a base model
mo = mo.genP(oir, inWei);  % generate the edge probability matrix
fprintf('%3.5fs\n',toc)

tic, fprintf('%-40s','Generating data ...')
mo = mo.genData;        % generate data (Adj. matrix "As" and the labels "c")
mo = mo.removeZeroDeg;  % remove zero degree nodes
fprintf('%3.5fs\n',toc)

% options for the init method and cpl/upl
init_opts = struct('verbose',false);
T = 20;
cpl_opts = struct('verbose',false,'delta_max',0.1, ...
                  'itr_num',T,'em_max',80,'track_err',false);

tic, fprintf('%-40s','Applying init. method (SCP) ...')  
% use spectral clustering with pertubation to initialize labels
[e init_dT] = initLabel5b(mo.As, mo.K, 'scp', init_opts);    
fprintf('%3.5fs\n',toc)
init_nmi = compErr(mo.c, e);


% apply cpl with init vector "e"
tic, fprintf('%-40s','Applying CPL ...') 
[chat, ~, cpl_dT, post] = ...
     cpl4c(mo.As, mo.K, e, mo.c, 'cpl', cpl_opts);
fprintf('%3.5fs\n',toc)
cpl_nmi = compErr(mo.c, chat);

% % estimate parameters (this is not going to be accurate for degree
% % corrected block model
% [pih Phat] = cplEstimParam('soft', mo.As, post, chat);
% [err_pri err_P] = compParamErr(mo.c, chat, mo.pri, pih, mo.P, Phat);

fprintf(1,'Init NMI = %3.2f\nCPL  NMI = %3.2f\n\n',init_nmi,cpl_nmi)
