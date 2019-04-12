addpath(fullfile('.','cpl_m_code'))

tic, fprintf('%-40s','Setting up the model ...')
n = 3000; % number of nodes
K = 6;    % number of communities
oir = 0.1;    % "O"ut-"I"n-"R"atio 
lambda = 22;    % average node degree

compErr = @(c,e) compMuI(compCM(c,e,K));    % use mutual info as a measure of error/sim.

inWei = ones(1,K);   % relative wieght of in-class probabilities
[lowVal lowProb] = deal(0.2,0.9); % An example of Degree corrected block model
%[lowVal lowProb] = deal(0.2,0); % An example of original block model


mo = dcBlkMod2(n,K,lambda, lowVal, lowProb); % create a base model
mo = mo.genP(oir, inWei);  % generate the edge probability matrix
fprintf('%3.5fs\n',toc)

tic, fprintf('%-40s','Generating data ...')
mo = mo.genData;        % generate data (Adj. matrix "As" and the labels "c")
mo = mo.removeZeroDeg;  % remove zero degree nodes
fprintf('%3.5fs\n',toc)

A = mo.As;
[chat, e] = find_labels(A,K);
% fprintf('%3.5fs\n',toc)
init_nmi = compErr(mo.c, e);
cpl_nmi = compErr(mo.c, chat);

% figure(1), clf, plot(mo.c,'.'), hold on, plot(e,'o')

% % estimate parameters (this is not going to be accurate for degree
% % corrected block model
% [pih Phat] = cplEstimParam('soft', mo.As, post, chat);
% [err_pri err_P] = compParamErr(mo.c, chat, mo.pri, pih, mo.P, Phat);

fprintf(1,'Init NMI = %3.2f\nCPL  NMI = %3.2f\n\n',init_nmi,cpl_nmi)

%%
Kvec = 2:10;
Klen = length(Kvec);
Tstat = zeros(Klen,2);
for r = 1:Klen
    Ktest = Kvec(r);
    Tstat(r,1) = PLTest(A, @find_labels, Ktest, 1, K);
    Tstat(r,2) = PLTest(A, @find_labels, Ktest, 2, K+1);
end

%%
figure(1), clf, hold on
plot(Kvec,Tstat(:,1),'ro-', 'MarkerFaceColor','r') %,'MarkerSize',10)
plot(Kvec,Tstat(:,2),'b.-','MarkerSize',10)

