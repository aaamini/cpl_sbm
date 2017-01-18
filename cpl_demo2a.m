tic, fprintf('%-40s','Setting up the model ...')
n = 80; % number of nodes
K = 2;    % number of communities
oir = 0.2;    % "O"ut-"I"n-"R"atio 
lambda = 7;    % average node degree

compErr = @(c,e) compMuI(compCM(c,e,K));    % use mutual info as a measure of error/sim.

inWei = [1 1];   % relative wieght of in-class probabilities
[lowVal lowProb] = deal(0.2,0); % An example of Degree corrected block model
%[lowVal lowProb] = deal(0.2,0); % An example of original block model


mo = dcBlkMod3(n,K,lambda, lowVal, lowProb); % create a base modl
mo = mo.genP(oir, inWei);  % generate the edge probability matrix
mo.P = [0.5  0.05;0.05 0.5];
fprintf('%3.5fs\n',toc)

tic, fprintf('%-40s','Generating data ...')
m = n/2;
mo = mo.genData([ones(1,m) 2*ones(1,m)]);        % generate data (Adj. matrix "As" and the labels "c")
%mo = mo.removeZeroDeg;  % remove pwidthzero degree nodes
fprintf('%3.5fs\n',toc)

%%
spy(mo.As)
print('-depsc','adj1.eps')

set(gca,'looseinset',[0 0 0 0])
mmwrite('presG1.mm',mo.As)
%%
% 
% % options for the init method and cpl/upl
% init_opts = struct('verbose',false);
% T = 20;
% cpl_opts = struct('verbose',false,'delta_max',0.1, ...
%                   'itr_num',T,'em_max',80,'track_err',false);
% 
% tic, fprintf('%-40s','Applying init. method (SCP) ...')  
% % use spectral clustering with pertubation to initialize labels
% [e init_dT] = initLabel5b(mo.As, mo.K, 'scp', init_opts);    
% fprintf('%3.5fs\n',toc)
% init_nmi = compErr(mo.c, e);
% 
% 
% % apply cpl with init vector "e"
% tic, fprintf('%-40s','Applying CPL ...') 
% [chat, ~, cpl_dT, post] = ...
%      cpl4c(mo.As, mo.K, e, mo.c, 'cpl', cpl_opts);
% fprintf('%3.5fs\n',toc)
% cpl_nmi = compErr(mo.c, chat);
% 
% % % estimate parameters (this is not going to be accurate for degree
% % % corrected block model
% % [pih Phat] = cplEstimParam('soft', mo.As, post, chat);
% % [err_pri err_P] = compParamErr(mo.c, chat, mo.pri, pih, mo.P, Phat);
% 
% fprintf(1,'Init NMI = %3.2f\nCPL  NMI = %3.2f\n\n',init_nmi,cpl_nmi)
% 
% mmwrite('presG1.mm',mo.As)
