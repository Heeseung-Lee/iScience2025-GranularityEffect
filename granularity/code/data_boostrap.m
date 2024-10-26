clear all; clc

addpath('./')
addpath('../data/')
cd('../result')

nIter       = 100;
nTB         = 10;
pCI         = 95;
lw          = 1.5;

%% trial number count

ipar        = load("data_ringpitch_seperate.mat");
nSubs([1 2])= size(ipar.pars,1);
data{1}     = ipar.pars(:,1);
data{2}     = ipar.pars(:,2);
ipar = load("data_ringpitch_alternate.mat");
data{3}     = ipar.pars;
nSubs(3)    = size(ipar.pars,1);
ipar = load("data_motor.mat");
data{4}     = ipar.pars;
nSubs(4)    = size(ipar.pars,1);
ipar = load("data_PreCur.mat");
data{5}     = ipar.pars;
nSubs(5)    = size(ipar.pars,1);
ipar = load("data_scaling.mat");
data{6}     = ipar.pars(12:end);
nSubs(6)    = size(ipar.pars(12:end),1);
data{7}     = ipar.pars(1:11);
nSubs(7)    = size(ipar.pars(1:11),1);

%%

features        = [1 2];
if isempty(dir(['results_RingPitch_nTB' num2str(nTB) '.mat'])) && isempty(dir(['results_RingPitch_nTB' num2str(nTB) '_compute.mat']))
    save(['results_RingPitch_nTB' num2str(nTB) '_compute.mat'],'features')
    
    pars = data(1:2);
    nSub = nSubs(1);

    r       = NaN(4*(2*nTB + 1),2);
    r_bst   = NaN(4*(2*nTB + 1),2,nIter);
    for ifeature = features
        X0s     = [];
        ids     = [];
        iss     = [];
        irts    = [];
        igs     = [];
        for icond = [1 2 3]
            
            d0  = [];
            s0  = [];
            rt0 = [];    
            for iSub = 1:nSub
                par         = pars{ifeature}{iSub};
                iIndR       = find((par.condition == icond) & (par.StairTrainTest==3));
                d0          = cat(2,d0,par.Chc(1:45,iIndR));
                s0          = cat(2,s0,par.Stm(1:45,iIndR));
                rt0         = cat(2,rt0,par.RT(1:45,iIndR));
            end

            d0          = zscore_HL(d0); % block-based z-score
            s0          = zscore_HL(s0);
            rt0         = zscore_HL(rt0);

            inR         = size(d0,2);
            id          = d0(:);
            is          = s0(:);
            irt         = rt0(:);
            for iTB = 1:nTB
                jd  = [NaN(iTB, inR); d0(1:end-iTB,:)];
                js  = [NaN(iTB, inR); s0(1:end-iTB,:)];
                jrt = [NaN(iTB, inR); rt0(1:end-iTB,:)];
                id  = cat(2,id,jd(:));
                is  = cat(2,is,js(:));
                irt = cat(2,irt,jrt(:));
            end
            jds         = [id is];
            ivalid      = ~isnan(sum(jds,2));
            id          = id(ivalid,:);
            is          = is(ivalid,:);
            irt         = irt(ivalid,:);

            ids         = cat(1,ids,id);
            iss         = cat(1,iss,is);
            irts        = cat(1,irts,irt);
            igs         = cat(1,igs,icond*ones(size(id,1),1));
        end
        ids             = zscore(ids); % final z-score
        iss             = zscore(iss);
        irts            = zscore(irts);
        X0              = [ids iss ids(:,2:end).*(igs-1) iss.*(igs-1) ids(:,2:end).*irts(:,2:end) iss.*irts ids(:,2:end).*irts(:,1) iss(:,2:end).*irts(:,1)];
        ir              = glmfit(X0(:,2:end),X0(:,1));
        r(:,ifeature)   = ir;

        inI     = size(X0,1);
        for iIter = 1:nIter
            fprintf('ifeature=%d/2 | iIter=%d/%d \n',ifeature,iIter,nIter)
            sInd                    = randi(inI,[inI,1]);
            ir                      = glmfit(X0(sInd,2:end),X0(sInd,1));
            r_bst(:,ifeature,iIter) = ir;
        end
    end

    save(['results_RingPitch_nTB' num2str(nTB) '.mat'],'r','r_bst')
    delete(['results_RingPitch_nTB' num2str(nTB) '_compute.mat'])    
end

%

dx = [-1.5 -0.5 0.5 1.5]*0.05;
colors = lines(4);
if ~isempty(dir(['results_RingPitch_nTB' num2str(nTB) '.mat']))
    load(['results_RingPitch_nTB' num2str(nTB) '.mat'])

    r_ci        = prctile(r_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    for ifeature = 1:2
        for ipar = 1:npar
            ibst    = squeeze(r_bst(ipar,ifeature,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,ifeature) = ip;
        end
        [~,~,~,pval_fdr(:,ifeature)] = fdr_bh(pval(:,ifeature));
    end
    
    % discard the bias term
    r           = r(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);
    
    figure(1)
    clf
    for ifeature = 1:2
        subplot(2,1,ifeature)
        hold on
        for icond = 4:-1:1
            if icond == 4
                iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB)-1);
                x       = [1:nTB nTB+2:(1+2*nTB)];
            else
                iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB));
                x       = (1:(1+2*nTB));
            end
            icoefs  = squeeze(r(iInd,ifeature));
            ici     = squeeze(r_ci(iInd,ifeature,:));
            ipval   = squeeze(pval_fdr(iInd,ifeature));
            for ixrange = 1:3
                switch ixrange
                    case 1
                        xInd = find(x<=nTB);
                    case 2
                        xInd = find(x==nTB+1);
                    case 3
                        xInd = find(x>nTB+1);
                end
                incoef  = length(xInd);
                for ix = 1:incoef
                    plot([x(xInd(ix)) x(xInd(ix))]+dx(icond),ici(xInd(ix),:),'_-','color',colors(icond,:))
                end
                plot(x(xInd)+dx(icond),icoefs(xInd),'o-','color',colors(icond,:),'MarkerFaceColor','w')     
                ip = ipval(xInd)<0.05;
                plot(x(xInd(ip))+dx(icond),icoefs(xInd(ip)),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))                                
            end
        end        
        plot([0 2*nTB+2],[0 0],'-k')
        xlim([0 2*nTB+2])
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1000 2000]/150)
    saveas(gcf,['Figure_1_RingPitch_' num2str(nTB) '.png'])
end

%%

if isempty(dir('results_GeneralMagnitude.mat')) && isempty(dir('results_GeneralMagnitude_compute.mat'))
    save('results_GeneralMagnitude_compute.mat','lw')

    pars = data{3};
    nSub = nSubs(3);

    coefs       = NaN(12,4); % 1:R2R, 2:P2R, 3:R2P, 4:P2P
    coefs_bst   = NaN(12,4,nIter);

    f0s     = [];
    f1s     = [];
    d0s     = [];
    s0s     = [];
    rt0s    = [];
    d1s     = [];
    s1s     = [];
    rt1s    = [];
    g1s     = [];
    for ig = [1 2] % {2AFC, 4AFC}
        for iSub = 1:nSub
            iIndR       = find((pars{iSub}.condition == ig) & (pars{iSub}.StairTrainTest==3));
            s0          = pars{iSub}.Stm(:,iIndR);
            jChc        = pars{iSub}.Chc(:,iIndR);
            rt0         = pars{iSub}.RT(:,iIndR);
            f0          = pars{iSub}.iFeature(:,iIndR); % 1:Ring, 2:Pitch
            nclass      = pars{iSub}.nclass(iIndR(1));
            inR         = size(s0,2);
            d0          = NaN(size(jChc));
            for ifeature = 1:2
                iInd        = (f0 == ifeature) & (jChc <= nclass*ifeature) & (jChc > nclass*(ifeature-1));
                d0(iInd)    = jChc(iInd) - nclass*(ifeature-1);
            end

            for if0 = 1:2
                iInd        = f0==if0;
                s0(iInd)    = zscore_HL(s0(iInd)); % block-based z-score
                d0(iInd)    = zscore_HL(d0(iInd));
                rt0(iInd)   = zscore_HL(rt0(iInd));
            end

            s1          = [NaN(1,inR); s0(1:end-1,:)];
            d1          = [NaN(1,inR); d0(1:end-1,:)];
            f1          = [NaN(1,inR); f0(1:end-1,:)];
            rt1         = [NaN(1,inR); rt0(1:end-1,:)];

            f0s     = cat(2,f0s,f0);
            f1s     = cat(2,f1s,f1);
            d0s     = cat(2,d0s,d0);
            s0s     = cat(2,s0s,s0);
            rt0s    = cat(2,rt0s,rt0);
            d1s     = cat(2,d1s,d1);
            s1s     = cat(2,s1s,s1);
            rt1s    = cat(2,rt1s,rt1);
            g1s     = cat(2,g1s,ig*ones(size(rt1)));
        end
    end

    cf = 1;
    for ifeature = 1:2
        for pfeature = 1:2
            jInd    = (f0s==ifeature) & (f1s==pfeature) & ~isnan(d0s + d1s);
            id0     = NaN(size(d0s));
            is0     = NaN(size(d0s));
            irt0    = NaN(size(d0s));
            id1     = NaN(size(d0s));
            is1     = NaN(size(d0s));
            irt1    = NaN(size(d0s));
            ig1     = NaN(size(d0s));

            id0(jInd)   = d0s(jInd);
            is0(jInd)   = s0s(jInd);
            irt0(jInd)  = rt0s(jInd);
            id1(jInd)   = d1s(jInd);
            is1(jInd)   = s1s(jInd);
            irt1(jInd)  = rt1s(jInd);
            ig1(jInd)   = g1s(jInd);

            ivalid      = ~isnan(sum([id0(:) is0(:) irt0(:) id1(:) is1(:) irt1(:)],2));

            id0s        = zscore(id0(ivalid)); % final z-score
            is0s        = zscore(is0(ivalid));
            irt0s       = zscore(irt0(ivalid));
            id1s        = zscore(id1(ivalid));
            is1s        = zscore(is1(ivalid));
            irt1s       = zscore(irt1(ivalid));
            ig1s        = ig1(ivalid);

            X0              = [id0s id1s is0s is1s id1s.*(ig1s-1) is0s.*(ig1s-1) is1s.*(ig1s-1) id1s.*irt0s is0s.*irt0s is1s.*irt0s id1s.*irt1s is1s.*irt1s];
            ir              = glmfit(X0(:,2:end),X0(:,1));
            coefs(:,cf)     = ir;

            inI         = size(X0,1);
            for iIter = 1:nIter
                fprintf('cfeature=%d/4 | iIter=%d/%d \n',cf,iIter,nIter)
                sInd                    = randi(inI,[inI,1]);
                ir                      = glmfit(X0(sInd,2:end),X0(sInd,1));
                coefs_bst(:,cf,iIter)   = ir;
            end

            cf = cf + 1;
        end
    end

    save('results_GeneralMagnitude.mat','coefs','coefs_bst')
    delete('results_GeneralMagnitude_compute.mat')
end

%

if ~isempty(dir('results_GeneralMagnitude.mat'))
    load('results_GeneralMagnitude.mat')

    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,4);
    pval_fdr    = NaN(npar,4);
    for iseq = 1:4
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [~,~,~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(2)
    
    clf
    colors = lines(4);
    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    for iseq = 1:4
        subplot(1,4,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
            if icond == 4
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end

    saveas(gcf,'Figure_2_GeneralMagnitude.png')

end

%%

if isempty(dir('results_Motor.mat')) && isempty(dir('results_Motor_compute.mat'))
    save('results_Motor_compute.mat','lw')

    pars = data{4};
    nSub = nSubs(4);

    coefs       = NaN(12,4);
    coefs_bst   = NaN(12,4,nIter);
    for ifeature = [1 2] % {ring, pitch}
        for ih = 1:2
            hs      = [];
            d0s     = [];
            s0s     = [];
            rt0s    = [];
            d1s     = [];
            s1s     = [];
            rt1s    = [];
            g1s     = [];
            for ig = [1 2] % {2AFC, 4AFC}
                for iSub = 1:nSub
                    par     = pars{iSub};
                    iIndR   = find((par.condition == ig) & (par.features == ifeature) & (par.StairTrainTest==3));
                    inR     = length(iIndR);
                    nclass  = par.nclass(iIndR(1));
                    s0      = par.Stm(:,iIndR);
                    d0      = par.Chc(:,iIndR);
                    rt0     = par.RT(:,iIndR);
                    s1      = [NaN(1,inR); s0(1:end-1,:)];
                    rt1     = [NaN(1,inR); rt0(1:end-1,:)];

                    % 2AFC: 1~2->left | 3~4->right
                    % 4AFC: 1~4->left | 5~8->right
                    iMotor                  = par.Motor(:,iIndR); % cued motor
                    iMiss                   = isnan(par.RT(:,iIndR));
                    HandUsed                = NaN(size(d0));
                    HandUsed(d0<=nclass)    = 1; % left
                    HandUsed(d0>nclass)     = 2; % right
                    HandMatched             = HandUsed == iMotor; % the hand-matched trials are only condisered

                    HandUsed(~HandMatched | iMiss)  = NaN;
                    pHandUsed                       = [NaN(1,inR); HandUsed(1:end-1,:)];
                    d0(~HandMatched | iMiss)        = NaN;
                    d1                              = [NaN(1,inR); d0(1:end-1,:)];
                    d0(HandUsed==2)                 = d0(HandUsed==2) - 2*ig;
                    d1(pHandUsed==2)                = d1(pHandUsed==2) - 2*ig;

                    iSameDiffHand                   = NaN(size(HandUsed));
                    iSameDiffHand(HandUsed == pHandUsed & ~isnan(pHandUsed) & ~isnan(HandUsed)) = 1;
                    iSameDiffHand(HandUsed ~= pHandUsed & ~isnan(pHandUsed) & ~isnan(HandUsed)) = 2;

                    hs      = cat(2,hs,iSameDiffHand);
                    d0s     = cat(2,d0s,zscore_HL(d0)); % block-based z-score
                    s0s     = cat(2,s0s,zscore_HL(s0));
                    rt0s    = cat(2,rt0s,zscore_HL(rt0));
                    d1s     = cat(2,d1s,zscore_HL(d1));
                    s1s     = cat(2,s1s,zscore_HL(s1));
                    rt1s    = cat(2,rt1s,zscore_HL(rt1));
                    g1s     = cat(2,g1s,ig*ones(size(rt1)));
                end
            end

            jInd    = hs==ih & ~isnan(d0s + d1s);

            d0s(~jInd)  = NaN;
            s0s(~jInd)  = NaN;
            rt0s(~jInd) = NaN;
            d1s(~jInd)  = NaN;
            s1s(~jInd)  = NaN;
            rt1s(~jInd) = NaN;

            ivalid      = ~isnan(sum([d0s(:) s0s(:) rt0s(:) d1s(:) s1s(:) rt1s(:)],2));

            d0s         = zscore(d0s(ivalid)); % final z-score
            s0s         = zscore(s0s(ivalid));
            rt0s        = zscore(rt0s(ivalid));
            d1s         = zscore(d1s(ivalid));
            s1s         = zscore(s1s(ivalid));
            rt1s        = zscore(rt1s(ivalid));
            g1s         = g1s(ivalid);

            X0                          = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt0s s0s.*rt0s s1s.*rt0s d1s.*rt1s s1s.*rt1s];
            ir                          = glmfit(X0(:,2:end),X0(:,1));
            coefs(:,2*(ifeature-1)+ih)  = ir;

            inI         = size(X0,1);
            for iIter = 1:nIter
                fprintf('ifeature=%d/2 | ih=%d/2 | iIter=%d/%d \n',ifeature,ih,iIter,nIter)
                sInd                                    = randi(inI,[inI,1]);
                ir                                      = glmfit(X0(sInd,2:end),X0(sInd,1));
                coefs_bst(:,2*(ifeature-1)+ih,iIter)    = ir;
            end
        end
    end

    save('results_Motor.mat','coefs','coefs_bst')
    delete('results_Motor_compute.mat')
end

%
if ~isempty(dir('results_Motor.mat'))
    load('results_Motor.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,4);
    h_fdr       = NaN(npar,4);
    pval_fdr    = NaN(npar,4);
    p_thre      = NaN(1,4);
    for iseq = 1:4
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(5)
    
    clf
    colors = lines(4);
    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    for iseq = 1:4
        subplot(1,4,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
            if icond == 4
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end

    saveas(gcf,'Figure_5_Motor.png')
end

%%

if isempty(dir('results_PreCur.mat')) && isempty(dir('results_PreCur_compute.mat'))
    save('results_PreCur_compute.mat','lw')

    pars = data{5};
    nSub = nSubs(5);

    coefs       = NaN(12,2);
    coefs_bst   = NaN(12,2,nIter);
    for icomb = 1:2
        d0s     = [];
        s0s     = [];
        rt0s    = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for iseq = (1:2)+2*(icomb-1) % 1:2(#)->2(cv), 2:8(#)->2(cv), 3:2(cv)->2(#), 4:2(cv)->8(#)
            for iSub = 1:nSub
                par     = pars{iSub};

                iIndR   = find((par.condition == (-mod(iseq,2)+2)) & (par.StairTrainTest==3));
                iRGind  = par.RGind;
                inR     = length(iIndR);

                id0     = par.Chc(:,iIndR);
                is0     = par.Stm(:,iIndR);
                irt0    = par.RT(:,iIndR);
                is1     = [NaN(1,inR); is0(1:end-1,:)];
                id1     = [NaN(1,inR); id0(1:end-1,:)];
                irt1    = [NaN(1,inR); irt0(1:end-1,:)];
                switch iseq
                    case 1 % 2(#)->2(cv)
                        iIndCurrent     = iRGind==1;
                    case 2 % 8(#)->2(cv)
                        iIndCurrent     = iRGind==1;
                    case 3 % 2(cv)->2(#)
                        iIndCurrent     = iRGind==2;
                    case 4 % 2(cv)->8(#)
                        iIndCurrent     = iRGind==2;
                end
                d0      = id0(iIndCurrent,:);
                s0      = is0(iIndCurrent,:);
                rt0     = irt0(iIndCurrent,:);
                d1      = id1(iIndCurrent,:);
                s1      = is1(iIndCurrent,:);
                rt1     = irt1(iIndCurrent,:);

                d0s     = cat(2,d0s,zscore_HL(d0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(s0));
                rt0s    = cat(2,rt0s,zscore_HL(rt0));
                d1s     = cat(2,d1s,zscore_HL(d1));
                s1s     = cat(2,s1s,zscore_HL(s1));
                rt1s    = cat(2,rt1s,zscore_HL(rt1));
                g1s     = cat(2,g1s,(2-mod(iseq,2))*ones(size(rt1)));       
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) rt0s(:) d1s(:) s1s(:) rt1s(:)],2));

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        rt0s        = zscore(rt0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        X0              = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt0s s0s.*rt0s s1s.*rt0s d1s.*rt1s s1s.*rt1s];
        ir              = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,icomb)  = ir;

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('PreCue=%d/2 | iIter=%d/%d \n',icomb,iIter,nIter)
            sInd                        = randi(inI,[inI,1]);
            ir                          = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,icomb,iIter)    = ir;
        end
    end

    save('results_PreCur.mat','coefs','coefs_bst')
    delete('results_PreCur_compute.mat')
end

%

if ~isempty(dir('results_PreCur.mat'))
    load('results_PreCur.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(6)
    clf
    colors = lines(4);
    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
            if icond == 4
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end
    saveas(gcf,'Figure_6_PreCur.png')
end

%%

if isempty(dir('results_Scaling.mat')) && isempty(dir('results_Scaling_compute.mat'))
    save('results_Scaling_compute.mat','lw')

    pars        = data{6};
    nSub        = nSubs(6);
    
    coefs           = NaN(9,2);
    coefs_bst       = NaN(9,2,nIter);
    for iScaleRepro = 1:2 % 1:scaling, 2:reproduction
        d0s     = [];
        s0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for isequence = 1:2 % 1:2->R, 2:8->R, 3:R->2, 4:R->8
            for iSub = 1:nSub
                par     = pars{iSub};

                iIndR   = find((par.condition == 2*(iScaleRepro-1) + (-mod(isequence,2)+2)) & (par.StairTrainTest==3));
                inR     = length(iIndR);
                iRGind  = par.RGind;
                
                d0      = par.Chc(:,iIndR);
                s0      = par.Stm(:,iIndR);
                rt0     = par.RT(:,iIndR);
                s1      = [NaN(1,inR); s0(1:end-1,:)];
                d1      = [NaN(1,inR); d0(1:end-1,:)];
                rt1     = [NaN(1,inR); rt0(1:end-1,:)];
                if iScaleRepro == 1
                    d0  = par.Scaling(:,iIndR);
                end

                switch isequence
                    case 1 % 2->R
                        iIndCurrent     = iRGind==2;
                    case 2 % 8->R
                        iIndCurrent     = iRGind==2;
                    case 3 % R->2
                        iIndCurrent     = iRGind==1;
                    case 4 % R->8
                        iIndCurrent     = iRGind==1;
                end
                id0     = d0(iIndCurrent,:);
                is0     = s0(iIndCurrent,:);
                id1     = d1(iIndCurrent,:);
                is1     = s1(iIndCurrent,:);
                irt1    = rt1(iIndCurrent,:);

                d0s     = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(is0));
                d1s     = cat(2,d1s,zscore_HL(id1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));            
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d1s(:) s1s(:) rt1s(:)],2));

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        X0                      = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        ir                      = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,iScaleRepro)    = ir;
        

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('iScaleRepro=%d/2 | iIter=%d/%d \n',iScaleRepro,iIter,nIter)
            sInd                            = randi(inI,[inI,1]);
            ir                              = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,iScaleRepro,iIter)  = ir;
        end
    end

    save('results_Scaling.mat','coefs','coefs_bst')
    delete('results_Scaling_compute.mat')
end

%

if ~isempty(dir('results_Scaling.mat'))
    load('results_Scaling.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(7)
    clf
    colors = lines(4);
    dx = [-1 0 1]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 3:-1:1)
            if icond == 3
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end
    saveas(gcf,'Figure_7_Scaling.png')
end

%%

if isempty(dir('results_Scaling_LongStimNoMask.mat')) && isempty(dir('results_Scaling_LongStimNoMask_compute.mat'))
    save('results_Scaling_LongStimNoMask_compute.mat','lw')

    pars = data{7};
    nSub = nSubs(7);

    coefs           = NaN(9,2);
    coefs_bst       = NaN(9,2,nIter);
    for iScaleRepro = 1:2 % 1:scaling, 2:reproduction
        d0s     = [];
        s0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for isequence = 1:2 % 1:2->R, 2:8->R, 3:R->2, 4:R->8
            for iSub = 1:nSub
                par     = pars{iSub};

                iIndR   = find((par.condition == 2*(iScaleRepro-1) + (-mod(isequence,2)+2)) & (par.StairTrainTest==3));
                inR     = length(iIndR);
                iRGind  = par.RGind;
                
                d0      = par.Chc(:,iIndR);
                s0      = par.Stm(:,iIndR);
                rt0     = par.RT(:,iIndR);
                s1      = [NaN(1,inR); s0(1:end-1,:)];
                d1      = [NaN(1,inR); d0(1:end-1,:)];
                rt1     = [NaN(1,inR); rt0(1:end-1,:)];
                if iScaleRepro == 1
                    d0  = par.Scaling(:,iIndR);
                end

                switch isequence
                    case 1 % 2->R
                        iIndCurrent     = iRGind==2;
                    case 2 % 8->R
                        iIndCurrent     = iRGind==2;
                    case 3 % R->2
                        iIndCurrent     = iRGind==1;
                    case 4 % R->8
                        iIndCurrent     = iRGind==1;
                end
                id0     = d0(iIndCurrent,:);
                is0     = s0(iIndCurrent,:);
                id1     = d1(iIndCurrent,:);
                is1     = s1(iIndCurrent,:);
                irt1    = rt1(iIndCurrent,:);

                d0s     = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(is0));
                d1s     = cat(2,d1s,zscore_HL(id1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));            
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d1s(:) s1s(:) rt1s(:)],2));

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        X0                      = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        ir                      = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,iScaleRepro)    = ir;        

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('iScaleRepro=%d/2 | iIter=%d/%d \n',iScaleRepro,iIter,nIter)
            sInd                            = randi(inI,[inI,1]);
            ir                              = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,iScaleRepro,iIter)  = ir;
        end
    end

    save('results_Scaling_LongStimNoMask.mat','coefs','coefs_bst')
    delete('results_Scaling_LongStimNoMask_compute.mat')
end

%

if ~isempty(dir('results_Scaling_LongStimNoMask.mat'))
    load('results_Scaling_LongStimNoMask.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(7)
    clf
    colors = lines(4);
    dx = [-1 0 1]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 3:-1:1)
            if icond == 3
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end

    saveas(gcf,'Figure_8_Scaling_LongStimNoMask.png')

end

%%

if isempty(dir('results_Scaling_ALL.mat')) && isempty(dir('results_Scaling_ALL_compute.mat'))
    save('results_Scaling_ALL_compute.mat','lw')

    [pars, nSub]    = RingGaborScalingReproductionDataLoad_ALL;

    coefs           = NaN(9,2);
    coefs_bst       = NaN(9,2,nIter);
    for iScaleRepro = 1:2 % 1:scaling, 2:reproduction
        d0s     = [];
        s0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for isequence = 1:2 % 1:2->R, 2:8->R, 3:R->2, 4:R->8
            for iSub = 1:nSub
                par     = pars{iSub};

                iIndR   = find((par.condition == 2*(iScaleRepro-1) + (-mod(isequence,2)+2)) & (par.StairTrainTest==3));
                inR     = length(iIndR);
                iRGind  = par.RGind;
                
                d0      = par.Chc(:,iIndR);
                s0      = par.Stm(:,iIndR);
                rt0     = par.RT(:,iIndR);
                s1      = [NaN(1,inR); s0(1:end-1,:)];
                d1      = [NaN(1,inR); d0(1:end-1,:)];
                rt1     = [NaN(1,inR); rt0(1:end-1,:)];
                if iScaleRepro == 1
                    d0  = par.Scaling(:,iIndR);
                end

                switch isequence
                    case 1 % 2->R
                        iIndCurrent     = iRGind==2;
                    case 2 % 8->R
                        iIndCurrent     = iRGind==2;
                    case 3 % R->2
                        iIndCurrent     = iRGind==1;
                    case 4 % R->8
                        iIndCurrent     = iRGind==1;
                end
                id0     = d0(iIndCurrent,:);
                is0     = s0(iIndCurrent,:);
                id1     = d1(iIndCurrent,:);
                is1     = s1(iIndCurrent,:);
                irt1    = rt1(iIndCurrent,:);

                d0s     = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(is0));
                d1s     = cat(2,d1s,zscore_HL(id1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));            
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d1s(:) s1s(:) rt1s(:)],2));

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        X0                      = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        ir                      = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,iScaleRepro)    = ir;        

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('iScaleRepro=%d/2 | iIter=%d/%d \n',iScaleRepro,iIter,nIter)
            sInd                            = randi(inI,[inI,1]);
            ir                              = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,iScaleRepro,iIter)  = ir;
        end
    end

    save('results_Scaling_ALL.mat','coefs','coefs_bst')
    delete('results_Scaling_ALL_compute.mat')
end

%

if ~isempty(dir('results_Scaling_ALL.mat'))
    load('results_Scaling_ALL.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(7)
    clf
    colors = lines(4);
    dx = [-1 0 1]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 3:-1:1)
            if icond == 3
                iInd    = (1+(icond-1)*3):(icond*3-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*3):(icond*3);
                x       = 1:3;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 3.5],[0 0],'-k')
        xlim([0.5 3.5])
    end

    saveas(gcf,'Figure_8_Scaling_ALL.png')

end

%%

if isempty(dir('results_Scaling_Reduce.mat')) && isempty(dir('results_Scaling_Reduce_compute.mat')) || 1
    save('results_Scaling_Reduce_compute.mat','lw')

    [pars, nSub]    = RingGaborScalingReproductionDataLoad;

    coefs           = NaN(6,2);
    coefs_bst       = NaN(6,2,nIter);
    for iReduce = 1:2 % 1:Original G8; 2:G8 is reduced to G2 
        d0s     = [];
        s0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for isequence = 1:2 % 1: G1=2; 2: G1=8 
            for iSub = 1:nSub
                par     = pars{iSub};
                if iReduce == 1
                    iG  = isequence;
                else
                    iG  = 2; % only G=8
                end
                iScaleRepro = 1; % fixed to scaling
                iIndR   = find((par.condition == 2*(iScaleRepro-1) + (-mod(iG,2)+2)) & (par.StairTrainTest==3));
                inR     = length(iIndR);
                iRGind  = par.RGind;
                
                d0      = par.Chc(:,iIndR);
                s0      = par.Stm(:,iIndR);
                rt0     = par.RT(:,iIndR);
                s1      = [NaN(1,inR); s0(1:end-1,:)];
                d1      = [NaN(1,inR); d0(1:end-1,:)];
                rt1     = [NaN(1,inR); rt0(1:end-1,:)];
                if iScaleRepro == 1
                    d0  = par.Scaling(:,iIndR);
                end

                iIndCurrent     = iRGind==2;
                id0     = d0(iIndCurrent,:);
                is0     = s0(iIndCurrent,:);
                id1     = d1(iIndCurrent,:);
                is1     = s1(iIndCurrent,:);
                irt1    = rt1(iIndCurrent,:);
                
                if isequence == 1 && iReduce == 2
                    id1(id1<5) = 1;
                    id1(id1>4) = 2;
                end

                d0s     = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(is0));
                d1s     = cat(2,d1s,zscore_HL(id1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));            
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d1s(:) s1s(:) rt1s(:)],2));
        
        %

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        % X0                      = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        X0                  = [d0s d1s s0s d1s.*(g1s-1) s0s.*(g1s-1) d1s.*rt1s];
        ir                  = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,iReduce)    = ir;                

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('iReduce=%d/2 | iIter=%d/%d \n',iReduce,iIter,nIter)
            sInd                        = randi(inI,[inI,1]);
            ir                          = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,iReduce,iIter)  = ir;
        end
    end
 
    save('results_Scaling_Reduce.mat','coefs','coefs_bst')
    delete('results_Scaling_Reduce_compute.mat')
end

%

if ~isempty(dir('results_Scaling_Reduce.mat'))
    load('results_Scaling_Reduce.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(7)
    clf
    colors = lines(4);
    dx = [-1 0 1]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 3:-1:1)
            if icond == 3
                iInd    = (1+(icond-1)*2):(icond*2-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*2):(icond*2);
                x       = 1:2;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 2.5],[0 0],'-k')
        xlim([0.5 2.5])
    end

    saveas(gcf,'Figure_8_Scaling_Reduce.png')

end


%

if isempty(dir('results_Scaling_LongStimNoMask_Reduce.mat')) && isempty(dir('results_Scaling_Reduce_LongStimNoMask_compute.mat')) || 1
    save('results_Scaling_LongStimNoMask_Reduce_compute.mat','lw')

    [pars, nSub]    = RingGaborScalingReproductionDataLoad_LongStimNoMask;

    coefs           = NaN(6,2);
    coefs_bst       = NaN(6,2,nIter);
    for iReduce = 1:2 % 1:Original G8; 2:G8 is reduced to G2 
        d0s     = [];
        s0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        for isequence = 1:2 % 1: G1=2; 2: G1=8 
            for iSub = 1:nSub
                par     = pars{iSub};
                if iReduce == 1
                    iG  = isequence;
                else
                    iG  = 2; % only G=8
                end
                iScaleRepro = 1; % fixed to scaling
                iIndR   = find((par.condition == 2*(iScaleRepro-1) + (-mod(iG,2)+2)) & (par.StairTrainTest==3));
                inR     = length(iIndR);
                iRGind  = par.RGind;
                
                d0      = par.Chc(:,iIndR);
                s0      = par.Stm(:,iIndR);
                rt0     = par.RT(:,iIndR);
                s1      = [NaN(1,inR); s0(1:end-1,:)];
                d1      = [NaN(1,inR); d0(1:end-1,:)];
                rt1     = [NaN(1,inR); rt0(1:end-1,:)];
                if iScaleRepro == 1
                    d0  = par.Scaling(:,iIndR);
                end

                iIndCurrent     = iRGind==2;
                id0     = d0(iIndCurrent,:);
                is0     = s0(iIndCurrent,:);
                id1     = d1(iIndCurrent,:);
                is1     = s1(iIndCurrent,:);
                irt1    = rt1(iIndCurrent,:);
                
                if isequence == 1 && iReduce == 2
                    id1(id1<5) = 1;
                    id1(id1>4) = 2;
                end

                d0s     = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                s0s     = cat(2,s0s,zscore_HL(is0));
                d1s     = cat(2,d1s,zscore_HL(id1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));            
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d1s(:) s1s(:) rt1s(:)],2));
        
        %

        d0s         = zscore(d0s(ivalid)); % final z-score
        s0s         = zscore(s0s(ivalid));
        d1s         = zscore(d1s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);

        % X0                      = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        X0                  = [d0s d1s s0s d1s.*(g1s-1) s0s.*(g1s-1) d1s.*rt1s];
        ir                  = glmfit(X0(:,2:end),X0(:,1));
        coefs(:,iReduce)    = ir;                

        inI         = size(X0,1);
        for iIter = 1:nIter
            fprintf('iReduce=%d/2 | iIter=%d/%d \n',iReduce,iIter,nIter)
            sInd                        = randi(inI,[inI,1]);
            ir                          = glmfit(X0(sInd,2:end),X0(sInd,1));
            coefs_bst(:,iReduce,iIter)  = ir;
        end
    end
 
    save('results_Scaling_LongStimNoMask_Reduce.mat','coefs','coefs_bst')
    delete('results_Scaling_LongStimNoMask_Reduce_compute.mat')
end

%

if ~isempty(dir('results_Scaling_LongStimNoMask_Reduce.mat'))
    load('results_Scaling_LongStimNoMask_Reduce.mat')
    
    r_ci        = prctile(coefs_bst,[0 100]+[1 -1]*(100-pCI)/2,3);
    npar        = size(r_ci,1);
    pval        = NaN(npar,2);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        for ipar = 1:npar
            ibst    = squeeze(coefs_bst(ipar,iseq,:));
            ip      = 2*min([sum(0<ibst) sum(0>ibst)])/nIter; % chance = 0
            pval(ipar,iseq) = ip;
        end
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    figure(7)
    clf
    colors = lines(4);
    dx = [-1 0 1]*0.05;
    for iseq = 1:2
        subplot(1,2,iseq)
        hold on        
        for icond = 2:-1:1 % RT was not presented (to show RT; 3:-1:1)
            if icond == 3
                iInd    = (1+(icond-1)*2):(icond*2-1);
                x       = [1 3];
            else
                iInd    = (1+(icond-1)*2):(icond*2);
                x       = 1:2;
            end
            icoefs  = squeeze(coefs(iInd,iseq));
            ici     = squeeze(r_ci(iInd,iseq,:));
            ipval   = squeeze(pval_fdr(iInd,iseq));
            for ixrange = 1:length(x)
                plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:))
            end
            plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w')
            ip = ipval<0.05;
            plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:))
        end
        plot([0.5 2.5],[0 0],'-k')
        xlim([0.5 2.5])
    end

    saveas(gcf,'Figure_8_Scaling_LongStimNoMask_Reduce.png')

end