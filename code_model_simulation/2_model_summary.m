clear all; clc; close all

dir0 = '/Users/heeseunglee/Documents/project/granularity/github';

addpath([dir0 '/subs'])

models = {'Standard','Adaptation','ThetaUpdate','ThetaConded','RhoConded'};
nmodel = length(models);

figure(1)
clf
for jmodel = 1:nmodel
    imodel      = models{jmodel};
    iBin        = 1; % 1: 4; 2: 6; 3: 8

    nT          = 1000;
    npar        = 10;

    sig_ms      = linspace(0.1,1,npar);
    sig_mms     = linspace(0.4,4,npar);
    sig_pri0s   = linspace(0.4,4,npar);
    sig_rs      = linspace(0.02,0.2,npar);
    dif_rs      = linspace(0.02,0.2,npar);
    gains       = linspace(0.1,1,npar);
    dif_ts      = linspace(0.2,2,npar);
    
    switch imodel
        case 'Standard'
            isavename   = sprintf([dir0 '/result/models/' imodel '/'...
                'ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_srs%.2f%.2f_drs%.2f%.2f_nT%d_npar%d'],...
                sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),sig_rs(1),sig_rs(end),dif_rs(1),dif_rs(end),nT,npar);
        case 'Adaptation'
            isavename   = sprintf([dir0 '/result/models/Adaptation/ms%.2f%.2f_gains%.2f%.2f_pri0%.2f%.2f_nT%d_npar%d'],...
                sig_ms(1),sig_ms(end),gains(1),gains(end),sig_pri0s(1),sig_pri0s(end),nT,npar);
        case 'ThetaUpdate'
            isavename   = sprintf([dir0 '/result/models/ThetaUpdate/ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_nT%d_npar%d'],...
                sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),nT,npar);
        case 'ThetaConded'
            isavename   = sprintf([dir0 '/result/models/ThetaConditioned/ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_ts%.2f%.2f_nT%d_npar%d'],...
                sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),dif_ts(1),dif_ts(end),nT,npar);
        case 'RhoConded'
            isavename   = sprintf([dir0 '/result/models/RhoConded/'...
                'ms%.2f%.2f_pri0%.2f%.2f_srs%.2f%.2f_drs%.2f%.2f_nT%d_npar%d'],...
                sig_ms(1),sig_ms(end),sig_pri0s(1),sig_pri0s(end),sig_rs(1),sig_rs(end),dif_rs(1),dif_rs(end),nT,npar);
    end
    cd(isavename)

    files       = dir([isavename '/Iter*.mat']);
    if isempty(dir([isavename '/Summary_' num2str(length(files)) '.mat']))
        acoefs = [];
        aperfs = [];
        for iIter = 1:length(files)
            fprintf('%d/%d\n',iIter,length(files))
            if iIter == 1
                load([isavename '/' files(iIter).name],'coefs','pars','perfs');
            else
                load([isavename '/' files(iIter).name],'coefs','perfs');
            end
            acoefs = cat(4,acoefs,coefs);
            aperfs = cat(2,aperfs,perfs);
        end
        acoefs      = mean(acoefs,4,'omitnan');
        param       = pars;
        perform     = mean(perfs,2,'omitnan');
        save([isavename '/Summary_' num2str(length(files)) '.mat'],'acoefs','param','perform')
    else
        load([isavename '/Summary_' num2str(length(files)) '.mat'])
    end

    load([dir0 '/data/results_Scaling_ALL_probit.mat'])
    PreCur = load([dir0 '/data/results_PreCur.mat']);
    coef_observe    = coefs; % bin = 4:10; 1st bin = 4
    ci_observe      = r_ci; 
    coef_observe    = [coef_observe(:,2) coef_observe(:,1)]; % to align with model [reproduction; scaling]
    ci_observe      = permute(ci_observe, [2, 1, 3]);
    ci_observe      = cat(1,ci_observe(2,:,:),ci_observe(1,:,:));
    perfs_observe   = mean(perfs,[2 3]);
    [~,~,ci]        = ttest(perfs_observe);
    iInd            = perform>ci(1) & perform<ci(2);
    acoefs          = acoefs(:,:,iInd); % [regressor, task, param]
    param_valid     = param(iInd,:);
    coef_valid      = acoefs;
    
    ierror          = squeeze(sum(abs(coef_observe([1 3 4 6],:) - acoefs([1 3 4 6],1:2,:)),[1 2]));
    coef_best       = acoefs(:,1:2,find(ierror==min(ierror),1))';
    param_best      = param_valid(find(ierror==min(ierror),1),:);
    
    coef_observe    = coef_observe';
    %% model simulation setting

    nT                  = 10000000;
    dr                  = 0.01;
    rrange0             = (-4+dr/2):dr:(4-dr/2);
    rrange              = (dr/2):dr:(1-dr/2);
    nrr                 = length(rrange);
    Gs                  = [2 8];
    nG                  = length(Gs);
    rhobounds           = cell(1,nG);
    for ig = 1:nG
        rhobounds{ig}   = linspace(0,1,Gs(ig)+1);
    end
    nStimBins   = [4 6 8];
    nBins       = length(nStimBins);
    dt          = 0.01;
    trange      = (-20+dt/2):dt:(20-dt/2);
    mu_pri0     = 0;

    %%

    subplot(nmodel,10,10*(jmodel-1)+1) % previous choice
    hold on
    xcoef = squeeze(coef_valid(1,1,:))';
    ycoef = squeeze(coef_valid(1,2,:))';
    plot(xcoef,ycoef,'.','color',[1 1 1]*0.6)
    plot([0 0],[-3 3],'k--','LineWidth',1)
    plot([-3 3],[0 0],'k--','LineWidth',1)
    x = coef_observe(1,1);
    y = coef_observe(2,1);
    plot(x,y,'or','MarkerFaceColor','r')
    plot(squeeze(ci_observe(1,1,:)),[y y],'|-r','LineWidth',1.4)
    plot([x,x],squeeze(ci_observe(2,1,:)),'_-r','LineWidth',1.4)
    x = coef_best(1,1);
    y = coef_best(2,1);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',8)
    % xlabel('reproduction \beta_c_1')
    % ylabel('scaling \beta_c_1')
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    set(gca,'xtick',[-0.12 0 0.12],'xticklabel',[],'ytick',[-0.12 0 0.12],'yticklabel',[])
    
    subplot(nmodel,10,10*(jmodel-1)+2) % granularity on previous choice
    hold on
    xcoef = squeeze(coef_valid(4,1,:))';
    ycoef = squeeze(coef_valid(4,2,:))';
    plot(xcoef,ycoef,'.','color',[1 1 1]*0.6)
    plot([0 0],[-3 3],'k--','LineWidth',1)
    plot([-3 3],[0 0],'k--','LineWidth',1)
    x = coef_observe(1,4);
    y = coef_observe(2,4);
    plot(x,y,'or','MarkerFaceColor','r')
    plot(squeeze(ci_observe(1,4,:)),[y y],'|-r','LineWidth',1.4)
    plot([x,x],squeeze(ci_observe(2,4,:)),'_-r','LineWidth',1.4)
    x = coef_best(1,4);
    y = coef_best(2,4);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',8)
    % xlabel('reproduction \beta_c_1_X_G')
    % ylabel('scaling \beta_c_1_X_G')
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    set(gca,'xtick',[-0.12 0 0.12],'xticklabel',[],'ytick',[-0.12 0 0.12],'yticklabel',[])

    subplot(nmodel,10,10*(jmodel-1)+3) % previous stimulus
    hold on
    xcoef = squeeze(coef_valid(3,1,:))';
    ycoef = squeeze(coef_valid(3,2,:))';
    plot(xcoef,ycoef,'.','color',[1 1 1]*0.6)
    plot([0 0],[-3 3],'k--','LineWidth',1)
    plot([-3 3],[0 0],'k--','LineWidth',1)
    x = coef_observe(1,3);
    y = coef_observe(2,3);
    plot(x,y,'or','MarkerFaceColor','r')
    plot(squeeze(ci_observe(1,3,:)),[y y],'|-r','LineWidth',1.4)
    plot([x,x],squeeze(ci_observe(2,3,:)),'_-r','LineWidth',1.4)
    x = coef_best(1,3);
    y = coef_best(2,3);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',8)
    % xlabel('reproduction \beta_s_1')
    % ylabel('scaling \beta_s_1')
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    set(gca,'xtick',[-0.12 0 0.12],'xticklabel',[],'ytick',[-0.12 0 0.12],'yticklabel',[])

    subplot(nmodel,10,10*(jmodel-1)+4) % granularity on previous stimulus
    hold on
    xcoef = squeeze(coef_valid(6,1,:))';
    ycoef = squeeze(coef_valid(6,2,:))';
    plot(xcoef,ycoef,'.','color',[1 1 1]*0.6)
    plot([0 0],[-3 3],'k--','LineWidth',1)
    plot([-3 3],[0 0],'k--','LineWidth',1)
    x = coef_observe(1,6);
    y = coef_observe(2,6);
    plot(x,y,'or','MarkerFaceColor','r')
    plot(squeeze(ci_observe(1,6,:)),[y y],'|-r','LineWidth',1.4)
    plot([x,x],squeeze(ci_observe(2,6,:)),'_-r','LineWidth',1.4)
    x = coef_best(1,6);
    y = coef_best(2,6);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',8)
    % ylabel('scaling \beta_s_1_X_G')
    % xlabel('reproduction \beta_s_1_X_G')
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    set(gca,'xtick',[-0.12 0 0.12],'xticklabel',[],'ytick',[-0.12 0 0.12],'yticklabel',[])
    
    if isempty(dir([isavename '/Summary_' num2str(length(files)) '_simulation_nT' num2str(nT) '.mat']))
        switch imodel
            case 'Standard'
                [t, indG, thetas, rhos, classes] = model_Standard(param_best, Gs, nT, nrr, rrange0, rrange, rhobounds, mu_pri0);
            case 'Adaptation'
                [t, indG, thetas, rhos, classes] = model_Adaptation(param_best, Gs, nT, trange, dt, rhobounds, mu_pri0);
            case 'ThetaUpdate'
                [t, indG, thetas, rhos, classes] = model_ThetaUpdate(param_best, Gs, nT, rhobounds, mu_pri0);
            case 'ThetaConded'
                [t, indG, thetas, rhos, classes] = model_ThetaConded(param_best, Gs, nT, nrr, trange, dt, rhobounds, mu_pri0);
            case 'RhoConded'
                [t, indG, thetas, rhos, classes] = model_RhoConded(param_best, Gs, nT, nrr, rrange0, rrange, rhobounds, mu_pri0);
        end
        [~, ~, ~, results] = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
        save([isavename '/Summary_' num2str(length(files)) '_simulation_nT' num2str(nT) '.mat'],'results')
    else
        load([isavename '/Summary_' num2str(length(files)) '_simulation_nT' num2str(nT) '.mat'])
    end
    
    %%
    gs      = results.g;
    nBin    = 4;
    xval    = linspace(0,1,nBin+1);
    fz      = 15;
    for imodeldata = 1:2
        for ivar = 1:3 % 1: reproduction; 2: scaling; 3: classification
            for ig = 1:2
                subplot(nmodel,12,12*(jmodel-1)+6+2*(ivar-1)+ig)
                hold on
                if imodeldata == 1
                    iInd        = gs==ig;
                    iyhat       = results.yhat(iInd,ivar);
                    iy          = results.y(iInd,ivar);
                    ix          = results.s1(iInd);
                    ic1         = results.d1(iInd);
                else
                    if ivar == 3
                        iInd    = PreCur.xy{1}.g==ig;
                        iyhat   = PreCur.xy{1}.yhat(iInd);
                        iy      = PreCur.xy{1}.y(iInd);
                        ix      = PreCur.xy{1}.s1(iInd);
                        ic1     = PreCur.xy{1}.d1(iInd);
                    else
                        iInd    = xy{3-ivar}.g==ig;
                        iyhat   = xy{3-ivar}.yhat(iInd);
                        iy      = xy{3-ivar}.y(iInd);
                        ix      = xy{3-ivar}.s1(iInd);
                        ic1     = xy{3-ivar}.d1(iInd);
                    end
                end

                [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
                [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
                ip0         = order/size(ix,1);
                jx          = NaN(size(ix));
                for iBin = 1:nBin
                    if iBin == nBin
                        iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                    else
                        iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                    end
                    jx(iInd) = iBin;
                end

                ing = length(unique(ic1));
                cl = jet(ing);
                for jc1 = 1:ing
                    ky      = NaN(nBin,1);
                    kyhat   = NaN(nBin,1);
                    for iBin = 1:nBin
                        iInd        = (jx==iBin) & (ic1==jc1);
                        jy          = iy(iInd);
                        if ~isempty(jy)
                            ky(iBin)    = mean(jy);
                            kyhat(iBin) = mean(iyhat(iInd));
                        end
                    end
                    ix = find(~isnan(ky));
                    if imodeldata == 1
                        plot(ix,kyhat(ix),'-','color',cl(jc1,:),'linewidth',1.2)
                    else
                        plot(ix,ky(ix),'o','color','k','markerfacecolor',cl(jc1,:))
                    end
                    if ivar == 1
                        set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[-1 0 1],'fontsize',fz)
                        ylim([-2 2])
                    elseif ivar == 2
                        set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                        ylim([0 1])
                    else
                        set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                        ylim([0 1])
                    end
                    xlim([0.5 nBin+0.5])
                end
            end
        end
    end
end

savedir = '/Users/heeseunglee/Documents/project/granularity/github/result/models/Summary';
if isempty(dir(savedir))
    mkdir(savedir)
end
cd(savedir)

set(gcf,'paperunit','inches','paperposition',[0 0 2400 1100]/150*1.2)
saveas(gcf,[savedir '/nT' num2str(nT) '.png'])