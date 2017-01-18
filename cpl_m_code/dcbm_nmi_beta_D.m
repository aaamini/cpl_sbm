warning off all

VERBOSE = false;
SAVE_FIG = true;
REMOVE_ZD = true;
PARAM_ESTIM = false; 

TAG = '_beta_gam00_w1510_E';
figName = {};
figName{1} = ['nmi' TAG];
figName{2} = ['pri' TAG];
figName{3} = ['P' TAG];
dataFile = ['nmi_par' TAG '.mat'];


%%
n = 1200; % number of nodes
K = 3; % number of communities

compErr = @(c,e) compMuI(compCM(c,e,K));    % use mutual info as a measure of error/sim.

oirV = linspace(0.0,0.2,20);   % "V"ector of "O"ut-"I"n-"R"atio 
threshV = zeros(size(oirV));

inWei = [1 5 10];
lambda = 7;
[lowVal lowProb] = deal(0.2,0); 
bmodel = dcBlkMod2(n,K,lambda, lowVal, lowProb); % create a base model 

Nrea = 1e3;
T = 20;
estim_type = 'soft'; % in case we are doing parameter estimateion 

init_methods =  {'DC', ...          
                 'SC', ...     
                 'SCP'};
            
numInitMtds = numel(init_methods);
             
            
init = struct('code', {'bi_deg','sc','scp'}, ...
              'opts', {...
                struct('degPert',0.01,'verbose',false), ...
                struct('verbose',false), ...
                struct('verbose',false) ...
                });
         
% CPL_FLAGS(k) : whether to start CPL/UPL with init_methods{k}
CPL_FLAGS = [true; false; true];

%%
numMtds = numInitMtds + sum(CPL_FLAGS*2);
methods = cell(numMtds,1);
k = 1;

for it = 1:numInitMtds
    
    methods{k} = init_methods{it};
    
    if CPL_FLAGS(it)
        k = k+1;
        methods{k} = ['CPL [' init_methods{it} ']'];
        k = k+1;
        methods{k} = ['UPL [' init_methods{it} ']'];
    end
    
    k = k + 1;
end

%%
[GdTV Gerr Gerr_pri Gerr_P] = deal( zeros(numMtds, numel(oirV)) );

cpl_opts = struct('verbose',false,'delta_max',0.1, ...
                  'itr_num',T,'em_max',80,'track_err',false);

for r = 1:numel(oirV)
    oir = oirV(r);
    
    fprintf(1,'\n OIR = %3.3f, gam = %2.2f,  w = (%s) [%3d of %3d  ]\n', ...
        oir, lowProb, num2str(inWei), r, numel(oirV)) 
    
    
    mob = bmodel;
    mob = mob.genP(oir, inWei);
      
    %threshV(r) = n * abs(P(1,1) - P(1,2)) / (K*sqrt(lambda));
    threshV(r) = min(mob.compThresh);
    
    [dTV err] = deal( zeros(numMtds, numel(oirV), Nrea) );        
    
    parfor t = 1:Nrea
        [tmp_err tmp_err_pri tmp_err_P tmp_dTV] = deal(zeros(numMtds,1));
    
        tic
        mo = mob;
        % generate block model
        mo = mo.genData;
%         mo = mo.perturbData(0.25);
        if REMOVE_ZD
            mo = mo.removeZeroDeg;
        end
        modGenT = toc;

        
        k = 1;
        for it = 1:numInitMtds
           [e tmp_dTV(k)] = initLabel5b(mo.As, mo.K, ...
                                       init(it).code, init(it).opts); 
           
           tmp_err(k) = compErr(mo.c, e);
           
           if PARAM_ESTIM
               [pih Phat] = cplEstimParam('hard', mo.As, nan(K,n), e);
               [tmp_err_pri(k) tmp_err_P(k)] = ...
                compParamErr(mo.c, e, mo.pri, pih, mo.P, Phat);
               
           end
           
           if CPL_FLAGS(it)
               for COND_UNCOND = ['c' 'u']
                   k = k+1;
                   [chat, ~, tmp_dTV(k), post] = ...
                       cpl4c(mo.As, mo.K, e, mo.c, COND_UNCOND, cpl_opts);
                   tmp_err(k) = compErr(mo.c, chat);

                   if PARAM_ESTIM
                       [pih Phat] = cplEstimParam(estim_type, mo.As, post, chat);
                       [tmp_err_pri(k) tmp_err_P(k)] = ...
                       compParamErr(mo.c, chat, mo.pri, pih, mo.P, Phat);
                   end
               end
               
%                k = k+1;
%                [chat, ~, tmp_dTV(k), post] = cpl2(mo.As, mo.K, e, mo.c, 'u', cpl_opts);
%                tmp_err(k) = compErr(mo.c, chat);
%                
%                if PARAM_ESTIM
%                    [pih Phat] = cplEstimParam(estim_type, mo.As, post, chat);
%                    [tmp_err_pri(k) tmp_err_P(k)] = ...
%                    compParamErr(mo.c, chat, mo.pri, pih, mo.P, Phat);
%                end
           end
           
           k = k+1;
           
        end % for-it
                
        
        %if (t == floor(Nrea/2))
        fprintf(1,'\n t = %3d  [OIR = %3.3f, gam = %2.2f,  w = (%s) [%3d of %3d  ] ]\n', ...
        t, oir, lowProb, num2str(inWei), r, numel(oirV)) 
            %fprintf(1,'%3d > \n',t)
            fprintf(1,'%10s|','ModGen')
            for k = 1:numMtds
                fprintf(1,'%10s|',methods{k})
            end
            fprintf(1,'\n')
            fprintf(1,'% 10.2f|',modGenT)
            for k = 1:numMtds
                fprintf(1,'% 10.2f|',tmp_dTV(k))
            end
            %fprintf(1,'\n  ModGen_dT = %2.2f\n',modGenT)
            
            %for k = 1:numMtds
            %    fprintf(1,['\n  ' methods{k} ' = %3.3f'], tmp_dTV(k))
            %end
            fprintf(1,'\n')
            fprintf(1,'  Total T = %2.2f\n',sum(tmp_dTV))
        %else
            % fprintf(1,'%3d ',t)
        %end
        
%         if (t == floor(Nrea/2))
%             fprintf(1,'%3d > ',t)
%             fprintf(1,'\n  ModGen_dT = %2.2f\n',modGenT)
%             
%             for k = 1:numMtds
%                 fprintf(1,['\n  ' methods{k} ' = %3.3f'], tmp_dTV(k))
%             end
%             fprintf(1,'\n')
%             fprintf(1,'  Total T = %2.2f\n',sum(tmp_dTV))
%         else
%             % fprintf(1,'%3d ',t)
%         end
        
        err(:,r,t) = tmp_err;
        dTV(:,r,t) = tmp_dTV;
        err_pri(:,r,t) = tmp_err_pri;
        err_P(:,r,t) = tmp_err_P;
    end
    
    for k = 1:numMtds
        Gerr(k,r) = mean(err(k,r,:));
        GdTV(k,r) = mean(dTV(k,r,:));        
        Gerr_pri(k,r) = mean(err_pri(k,r,:));
        Gerr_P(k,r) = mean(err_P(k,r,:));
    end
    
end

%%
save(dataFile)

%%
wStr = ['w = (' sprintf('%d, ', inWei(1:end-1)) ...
                           sprintf('%d)',inWei(end))];
textStr = sprintf(...
          ['%s\n\\gamma = %1.1f\n\\lambda = %1.0f\nn = %4d\n' ...
          'K = %2d\nN_{rea} = %d\n'], ...
          wStr, lowProb, lambda, n, K, Nrea);

% textStr = sprintf(...
%           ['\\lambda = %1.0f,  \\gamma = %1.1f,  n = %4d,  ' ...
%           'K = %2d,  N_{rep} = %d\n'], ...
%           lambda, lowProb, n, K, Nrea);
%       
% 
%       
% textStr = [textStr 'w = (' sprintf('%d, ', inWei(1:end-1)) ...
%                            sprintf('%d)',inWei(end))];
      
cmap = colormap(jet);

% style = {'-','-','-', ...
%          's-', ...
%          '*--','*--','*--'};

% color = {'k', 'r', 'b', ...
%          'g', ...
%          [0.2 1 0.2], 'r', 'b'};
C1 = [0.2 1 0.2];
C2 = [0 .5 .6];
C3 = [.1 .3 .8];

style = {'-.', '.-', '-', '^-', '*--','s--','o--'};
color = {'k' , 'r', C2,   C1,   'g',  'r',  C3};  
mSize = {   5,   5,   5,    5,    5 ,    8,    5}; % marker size
     
IDX = [1:7]; 
numIDX = numel(IDX);

yLABs = {'Normalized mutual information', ...
          '2-norm error for \pi', ...
          'log_{10} || \Delta P||_2'};
GGERR = {Gerr, Gerr_pri, log10(Gerr_P)};
yLIMs = {[0 1],[0 0.8],[-2.5 -0.3]};
legLOC = {'NorthEast','NorthWest','NorthEast'};
txtLOC = {[0.05,0.48],[0.3,0.26],[0.05,0.96]};

if PARAM_ESTIM
    nFig = 3;
else
    nFig = 1;
end

close all
for FIG = 1:nFig
    set(figure(FIG), 'position', [300 380 630 530]), clf, 
    ax1 = axes('position',[0.1 0.1 0.85 0.8]); hold on

    ph = zeros(numIDX,1);
    for id = 1:numIDX
        k = IDX(id);
        ph(id) = plot(oirV, GGERR{FIG}(k,:), style{id}, ...
            'color', color{id},'MarkerSize',mSize{id});
    end

    set(ph,'LineWidth',2)
    hYLabel = ylabel(yLABs{FIG});
    hXLabel = xlabel('\beta');
    set([hXLabel, hYLabel], 'FontName', 'Arial', 'FontSize',14, ...
        'FontWeight','bold')

    axis([[min(oirV), max(oirV)] yLIMs{FIG}])

    th = text(txtLOC{FIG}(1),txtLOC{FIG}(2),textStr,'FontSize',15, ...
        'units','normalized','FontWeight','bold');


    set(th,'VerticalAlignment','top')
    hLeg = legend(ph, methods(IDX),...
                      'location',legLOC{FIG},'FontSize',13,'FontWeight','bold');
%    set(hLeg,'Position',[0.72 0.21 0.22 0.2])  
    set(hLeg,'Box','off')
    set(gca,'FontSize',14,'FontWeight','bold')
    
    if SAVE_FIG
        print('-depsc2', [figName{FIG} '.eps'])
    end
end
warning on all
