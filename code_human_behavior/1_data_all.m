clear all; clc

nTB         = 10;
nTB_indiv   = 7;

downloadfolder = '/Users/heeseunglee/Documents/project/granularity/github';
savedir     = [downloadfolder '/result/humans'];
if isempty(dir(savedir))
    mkdir(savedir)
end
cd(savedir)

addpath([downloadfolder '/subs'])
addpath([downloadfolder '/data'])

ipar        = load("data_ringpitch_seperate.mat");
nSub([1 2]) = size(ipar.pars,1);
pars{1}     = ipar.pars(:,1);
pars{2}     = ipar.pars(:,2);
ipar        = load("data_ringpitch_alternate.mat");
pars{3}     = ipar.pars;
nSub(3)     = size(ipar.pars,1);
ipar        = load("data_motor.mat");
pars{4}     = ipar.pars;
nSub(4)     = size(ipar.pars,1);
ipar        = load("data_PreCur.mat");
pars{5}     = ipar.pars;
nSub(5)     = size(ipar.pars,1);
ipar        = load("data_scaling.mat");
pars{6}     = ipar.pars(12:end);
nSub(6)     = size(ipar.pars(12:end),1);
pars{7}     = ipar.pars(1:11);
nSub(7)     = size(ipar.pars(1:11),1);
pars{8}     = ipar.pars;
nSub(8)     = size(ipar.pars,1);
opts        = detectImportOptions('BehavData.csv');
opts.SelectedVariableNames = {'task','subj','testrun_total','iti','granule','Stm','Chc','RT'};
idat        = readtable('BehavData.csv',opts);
iInd        = idat.task==3 & ~isnan(idat.testrun_total);
data{9}     = idat(iInd,:);
nSubs(9)    = 16;

%%
for z = 1
    %% RingPitch average

    features        = [1 2];

    bis     = [1 2 4];
    r       = [];
    pval    = [];
    r_ci    = [];
    xy = cell(1,2);
    for ifeature = features
        X0s     = [];
        ids     = [];
        iss     = [];
        irts    = [];
        igs     = [];
        iys     = [];
        id1s    = [];
        for icond = [1 2 3]

            d0  = [];
            s0  = [];
            rt0 = [];
            for iSub = 1:nSub(ifeature)
                par         = pars{ifeature}{iSub};
                iIndR       = find((par.condition == icond) & (par.StairTrainTest==3));
                d0          = cat(2,d0,par.Chc(1:45,iIndR));
                s0          = cat(2,s0,par.Stm(1:45,iIndR));
                rt0         = cat(2,rt0,par.RT(1:45,iIndR));
            end

            y0          = NaN(size(d0));
            y0(d0>bis(icond))   = 1;
            y0(d0<=bis(icond))  = 0;
            y0          = y0(:);

            inR         = size(d0,2);
            d1          = [NaN(1,inR); d0(1:end-1,:)];
            d1          = d1(:);

            [id01,id02,id03] = binarizer(d0,ones(size(d0))*bis(icond)*2);

            id01        = zscore_HL(id01); % block-based z-score
            id02        = zscore_HL(id02);
            id03        = zscore_HL(id03);
            s0          = zscore_HL(s0);
            rt0         = zscore_HL(rt0);

            kd01        = [];
            kd02        = [];
            kd03        = [];
            is          = s0(:);
            irt         = rt0(:);
            for iTB = 1:nTB
                jd01    = [NaN(iTB, inR); id01(1:end-iTB,:)];
                jd02    = [NaN(iTB, inR); id02(1:end-iTB,:)];
                jd03    = [NaN(iTB, inR); id03(1:end-iTB,:)];
                js      = [NaN(iTB, inR); s0(1:end-iTB,:)];
                jrt     = [NaN(iTB, inR); rt0(1:end-iTB,:)];
                kd01    = cat(2,kd01,jd01(:));
                kd02    = cat(2,kd02,jd02(:));
                kd03    = cat(2,kd03,jd03(:));
                is      = cat(2,is,js(:));
                irt     = cat(2,irt,jrt(:));
            end
            if icond == 1
                kd02    = zeros(size(kd01));
                kd03    = zeros(size(kd01));
            elseif icond == 2
                kd03    = zeros(size(kd01));
            end
            jds         = [kd01 kd02 kd03];
            ivalid      = ~isnan(sum([jds is irt],2));
            id          = jds(ivalid,:);
            is          = is(ivalid,:);
            irt         = irt(ivalid,:);
            iy          = y0(ivalid,:);
            id1         = d1(ivalid,:);

            ids         = cat(1,ids,id);
            iss         = cat(1,iss,is);
            irts        = cat(1,irts,irt);
            iys         = cat(1,iys,iy);
            id1s        = cat(1,id1s,id1);
            igs         = cat(1,igs,icond*ones(size(id,1),1));
        end
        iInd            = ids==0;
        kds             = NaN(size(ids));
        kds(~iInd)      = ids(~iInd);
        kds             = zscore_HL(kds); % final z-score
        kds(iInd)       = 0;
        ids             = kds;
        iss             = zscore(iss);
        irts            = zscore(irts);
        
        ids_orig        = ids(:,1:nTB);
        ids_cont        = ids(:,nTB+1:end);
        X0_orig         = [iys ids_orig iss ids_orig.*(igs-1) iss.*(igs-1) ids_orig.*irts(:,2:end) iss.*irts ids_orig.*irts(:,1) iss(:,2:end).*irts(:,1)];
        X0_cont         = [ids_cont ids_cont(:,1:nTB).*(igs-1) ids_cont.*repmat(irts(:,2:end),1,2) ids_cont.*irts(:,1)];
        X0              = [X0_orig X0_cont];
        [ir,~,istat]        = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
        r(:,ifeature)       = ir;
        pval(:,ifeature)    = istat.p;
        r_ci(:,ifeature,:)  = istat.se*1.96*[-1 1]+ir;

        xy{ifeature}.y      = iys;
        xy{ifeature}.s0     = iss(:,1);
        xy{ifeature}.s1     = iss(:,2);
        xy{ifeature}.d1     = id1s;
        xy{ifeature}.g      = igs;
        xy{ifeature}.yhat   = glmval(ir,X0(:,2:end),'probit');
    end
    
    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    colors = lines(4);
    
    npar        = size(r_ci,1);
    pval_fdr    = NaN(npar,2);
    for ifeature = 1:2
        [~,~,~,pval_fdr(:,ifeature)] = fdr_bh(pval(:,ifeature));
    end

    % discard the bias term
    r           = r(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    coefs = r;
    save(['results_RingPitch_nTB' num2str(nTB) '.mat'],'coefs','r_ci','pval_fdr','xy','pval');

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

    %% RingPitch individual

    features        = [1 2];

    bis     = [1 2 4];
    r       = [];
    irs     = [];
    pval    = [];
    r_ci    = [];
    xy      = cell(nSub(1),2);
    for ifeature = features
        X0s     = [];
        ids     = [];
        iss     = [];
        irts    = [];
        igs     = [];
        iys     = [];
        id1s    = [];
        sjs     = [];
        rs      = [];
        for icond = [1 2 3]

            d0  = [];
            s0  = [];
            rt0 = [];
            sj  = [];
            for iSub = 1:nSub(ifeature)
                par         = pars{ifeature}{iSub};
                iIndR       = find((par.condition == icond) & (par.StairTrainTest==3));
                d0          = cat(2,d0,par.Chc(1:45,iIndR));
                s0          = cat(2,s0,par.Stm(1:45,iIndR));
                rt0         = cat(2,rt0,par.RT(1:45,iIndR));
                sj          = cat(2,sj,iSub*ones(size(par.RT(1:45,iIndR))));
            end

            y0          = NaN(size(d0));
            y0(d0>bis(icond))   = 1;
            y0(d0<=bis(icond))  = 0;
            y0          = y0(:);

            inR         = size(d0,2);
            d1          = [NaN(1,inR); d0(1:end-1,:)];
            d1          = d1(:);

            [id01,id02,id03] = binarizer(d0,ones(size(d0))*bis(icond)*2);

            id01        = zscore_HL(id01); % block-based z-score
            id02        = zscore_HL(id02);
            id03        = zscore_HL(id03);
            s0          = zscore_HL(s0);
            rt0         = zscore_HL(rt0);

            inR         = size(d0,2);
            kd01        = [];
            kd02        = [];
            kd03        = [];
            is          = s0(:);
            irt         = rt0(:);
            sj          = sj(:);
            for iTB = 1:nTB_indiv
                jd01    = [NaN(iTB, inR); id01(1:end-iTB,:)];
                jd02    = [NaN(iTB, inR); id02(1:end-iTB,:)];
                jd03    = [NaN(iTB, inR); id03(1:end-iTB,:)];
                js      = [NaN(iTB, inR); s0(1:end-iTB,:)];
                jrt     = [NaN(iTB, inR); rt0(1:end-iTB,:)];
                kd01    = cat(2,kd01,jd01(:));
                kd02    = cat(2,kd02,jd02(:));
                kd03    = cat(2,kd03,jd03(:));
                is      = cat(2,is,js(:));
                irt     = cat(2,irt,jrt(:));
            end
            if icond == 1
                kd02    = zeros(size(kd01));
                kd03    = zeros(size(kd01));
            elseif icond == 2
                kd03    = zeros(size(kd01));
            end
            jds         = [kd01 kd02 kd03];
            ivalid      = ~isnan(sum([jds is irt],2));
            id          = jds(ivalid,:);
            is          = is(ivalid,:);
            irt         = irt(ivalid,:);
            iy          = y0(ivalid,:);
            id1         = d1(ivalid,:);
            sj          = sj(ivalid,:);

            ids         = cat(1,ids,id);
            iss         = cat(1,iss,is);
            irts        = cat(1,irts,irt);
            iys         = cat(1,iys,iy);
            id1s        = cat(1,id1s,id1);
            igs         = cat(1,igs,icond*ones(size(id,1),1));
            sjs         = cat(1,sjs,sj);
        end
        for iSub = 1:nSub(1)
            sInd            = sjs==iSub;
            sds             = ids(sInd,:);
            iInd            = sds==0;
            kds             = NaN(size(sds));
            kds(~iInd)      = sds(~iInd);
            kds             = zscore_HL(kds); % final z-score
            kds(iInd)       = 0;
            jds             = kds;
            jss             = zscore(iss(sInd,:));
            jrts            = zscore(irts(sInd,:));

            ids_orig        = jds(:,1:nTB_indiv);
            ids_cont        = jds(:,nTB_indiv+1:end);

            jy              = iys(sInd);
            jd1             = id1s(sInd);
            jgs             = igs(sInd);

            X0_orig         = [jy ids_orig jss ids_orig.*(jgs-1) jss.*(jgs-1) ids_orig.*jrts(:,2:end) jss.*jrts ids_orig.*jrts(:,1) jss(:,2:end).*jrts(:,1)];
            X0_cont         = [ids_cont ids_cont(:,1:nTB_indiv).*(jgs-1) ids_cont.*repmat(jrts(:,2:end),1,2) ids_cont.*jrts(:,1)];
            X0              = [X0_orig X0_cont];
            ir              = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
            rs(:,iSub)      = ir;

            xy{iSub,ifeature}.y     = jy;
            xy{iSub,ifeature}.s0    = jss(:,1);
            xy{iSub,ifeature}.s1    = jss(:,2);
            xy{iSub,ifeature}.d1    = jd1;
            xy{iSub,ifeature}.g     = jgs;
            xy{iSub,ifeature}.yhat  = glmval(ir,X0(:,2:end),'probit');
        end
        rs = rs';
        [~,ip,ci] = ttest(rs);

        r(:,ifeature)       = mean(rs);
        irs(:,ifeature,:)   = rs';
        pval(:,ifeature)    = ip;
        r_ci(:,ifeature,:)  = ci';
    end

    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    colors = lines(4);

    npar        = size(r_ci,1);
    pval_fdr    = NaN(npar,2);
    for ifeature = 1:2
        [~,~,~,pval_fdr(:,ifeature)] = fdr_bh(pval(:,ifeature));
    end

    % discard the bias term
    r           = r(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);
    irs         = irs(2:end,:,:);

    coefs = r;
    save(['results_RingPitch_nTB' num2str(nTB_indiv) '_Individual.mat'],'coefs','r_ci','pval_fdr','xy','pval','irs');

    figure(11)
    clf
    for ifeature = 1:2
        subplot(2,1,ifeature)
        hold on
        for icond = 4:-1:1
            if icond == 4
                iInd    = (1+(icond-1)*(1+2*nTB_indiv)):(icond*(1+2*nTB_indiv)-1);
                x       = [1:nTB_indiv nTB_indiv+2:(1+2*nTB_indiv)];
            else
                iInd    = (1+(icond-1)*(1+2*nTB_indiv)):(icond*(1+2*nTB_indiv));
                x       = (1:(1+2*nTB_indiv));
            end
            icoefs  = squeeze(r(iInd,ifeature));
            ici     = squeeze(r_ci(iInd,ifeature,:));
            ipval   = squeeze(pval_fdr(iInd,ifeature));
            for ixrange = 1:3
                switch ixrange
                    case 1
                        xInd = find(x<=nTB_indiv);
                    case 2
                        xInd = find(x==nTB_indiv+1);
                    case 3
                        xInd = find(x>nTB_indiv+1);
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
        plot([0 2*nTB_indiv+2],[0 0],'-k')
        xlim([0 2*nTB_indiv+2])
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1000 2000]/150)
    saveas(gcf,['Figure_1_RingPitch_' num2str(nTB_indiv) '_Individual.png'])

    %% GeneralMagnitude average

    coefs   = []; % 1:R2R, 2:P2R, 3:R2P, 4:P2P
    pval    = [];
    r_ci    = [];
    xy      = cell(2,2);

    bis     = [1 2];
    f0s     = [];
    f1s     = [];
    d0s     = [];
    s0s     = [];
    rt0s    = [];
    d01s    = [];
    d02s    = [];
    s1s     = [];
    rt1s    = [];
    g1s     = [];
    jd1s    = [];
    for ig = [1 2] % {2AFC, 4AFC}
        for iSub = 1:nSub(3)
            iIndR       = find((pars{3}{iSub}.condition == ig) & (pars{3}{iSub}.StairTrainTest==3));
            s0          = pars{3}{iSub}.Stm(:,iIndR);
            jChc        = pars{3}{iSub}.Chc(:,iIndR);
            rt0         = pars{3}{iSub}.RT(:,iIndR);
            f0          = pars{3}{iSub}.iFeature(:,iIndR); % 1:Ring, 2:Pitch
            nclass      = pars{3}{iSub}.nclass(iIndR(1));
            inR         = size(s0,2);
            d0          = NaN(size(jChc));
            for ifeature = 1:2
                iInd        = (f0 == ifeature) & (jChc <= nclass*ifeature) & (jChc > nclass*(ifeature-1));
                d0(iInd)    = jChc(iInd) - nclass*(ifeature-1);
            end
            y0              = NaN(size(d0));
            y0(d0>bis(ig))  = 1;
            y0(d0<=bis(ig)) = 0;

            [id01,id02,id03] = binarizer(d0,ones(size(d0))*bis(ig)*2);

            inR         = size(d0,2);
            id1         = [NaN(1,inR); d0(1:end-1,:)];

            for if0 = 1:2
                iInd        = f0==if0;
                imat        = NaN(size(s0));

                imat(iInd)  = s0(iInd);
                imat        = zscore_HL(imat); % block-based z-score
                s0(iInd)    = imat(iInd);

                imat(iInd)  = id01(iInd);
                imat        = zscore_HL(imat); % block-based z-score
                id01(iInd)  = imat(iInd);

                if ig == 2
                    imat(iInd)  = id02(iInd);
                    imat        = zscore_HL(imat); % block-based z-score
                    id02(iInd)  = imat(iInd);
                else
                    id02(iInd)  = 0;
                end

                imat(iInd)  = rt0(iInd);
                imat        = zscore_HL(imat); % block-based z-score
                rt0(iInd)   = imat(iInd);
            end

            s1          = [NaN(1,inR); s0(1:end-1,:)];
            jd01        = [NaN(1, inR); id01(1:end-1,:)];
            jd02        = [NaN(1, inR); id02(1:end-1,:)];
            f1          = [NaN(1,inR); f0(1:end-1,:)];
            rt1         = [NaN(1,inR); rt0(1:end-1,:)];

            f0s     = cat(2,f0s,f0);
            f1s     = cat(2,f1s,f1);
            d0s     = cat(2,d0s,y0);
            s0s     = cat(2,s0s,s0);
            rt0s    = cat(2,rt0s,rt0);
            d01s    = cat(2,d01s,jd01);
            d02s    = cat(2,d02s,jd02);
            s1s     = cat(2,s1s,s1);
            rt1s    = cat(2,rt1s,rt1);
            g1s     = cat(2,g1s,ig*ones(size(rt1)));
            jd1s    = cat(2,jd1s,id1);
        end
    end

    cf = 1;
    X0s = cell(2,2);
    for ifeature = 1:2
        for pfeature = 1:2
            jInd    = (f0s==ifeature) & (f1s==pfeature) & ~isnan(d0s + d01s + d02s);
            id0     = NaN(size(d0s));
            is0     = NaN(size(d0s));
            irt0    = NaN(size(d0s));
            id01    = NaN(size(d0s));
            id02    = NaN(size(d0s));
            is1     = NaN(size(d0s));
            irt1    = NaN(size(d0s));
            ig1     = NaN(size(d0s));
            jd1     = NaN(size(d0s));

            id0(jInd)   = d0s(jInd);
            is0(jInd)   = s0s(jInd);
            irt0(jInd)  = rt0s(jInd);
            id01(jInd)  = d01s(jInd);
            id02(jInd)  = d02s(jInd);
            is1(jInd)   = s1s(jInd);
            irt1(jInd)  = rt1s(jInd);
            ig1(jInd)   = g1s(jInd);
            jd1(jInd)   = jd1s(jInd);

            ivalid      = ~isnan(sum([id0(:) is0(:) irt0(:) id01(:) id02(:) is1(:) irt1(:)],2));

            id0s        = id0(ivalid);
            is0s        = zscore(is0(ivalid)); % final z-score
            irt0s       = zscore(irt0(ivalid));
            id01s       = zscore(id01(ivalid));
            uds         = id02(ivalid);
            iInd        = uds == 0;
            kds         = NaN(size(uds));
            kds(~iInd)  = uds(~iInd);
            kds         = zscore_HL(kds);
            kds(iInd)   = 0;
            id02s       = kds;
            is1s        = zscore(is1(ivalid));
            irt1s       = zscore(irt1(ivalid));
            ig1s        = ig1(ivalid);
            kd1s        = jd1(ivalid);

            ids_orig    = id01s;
            ids_cont    = id02s;

            X0_orig         = [id0s ids_orig is0s is1s ids_orig.*(ig1s-1) is0s.*(ig1s-1) is1s.*(ig1s-1) ids_orig.*irt0s is0s.*irt0s is1s.*irt0s ids_orig.*irt1s is1s.*irt1s];
            X0_cont         = [ids_cont ids_cont.*irt0s ids_cont.*irt1s];
            X0              = [X0_orig X0_cont];
            X0s{ifeature,pfeature} = X0;
            [ir,id,istat]   = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
            coefs(:,cf)     = ir;
            pval(:,cf)      = istat.p;
            r_ci(:,cf,:)    = istat.se*1.96*[-1 1]+ir;

            xy{ifeature,pfeature}.y      = id0s;
            xy{ifeature,pfeature}.s0     = is0s;
            xy{ifeature,pfeature}.s1     = is1s;
            xy{ifeature,pfeature}.d1     = kd1s;
            xy{ifeature,pfeature}.g      = ig1s;
            xy{ifeature,pfeature}.yhat   = glmval(ir,X0(:,2:end),'probit');
            cf = cf + 1;
        end
    end

    for samediff = 1:2
        if samediff == 1
            X0 = cat(1,X0s{1,1},X0s{2,2});
            id0s = cat(1,xy{1,1}.y,xy{2,2}.y);
            is0s = cat(1,xy{1,1}.s0,xy{2,2}.s0);
            is1s = cat(1,xy{1,1}.s1,xy{2,2}.s1);
            kd1s = cat(1,xy{1,1}.d1,xy{2,2}.d1);
            ig1s = cat(1,xy{1,1}.g,xy{2,2}.g);
        else
            X0 = cat(1,X0s{1,2},X0s{2,1});
            id0s = cat(1,xy{1,2}.y,xy{2,1}.y);
            is0s = cat(1,xy{1,2}.s0,xy{2,1}.s0);
            is1s = cat(1,xy{1,2}.s1,xy{2,1}.s1);
            kd1s = cat(1,xy{1,2}.d1,xy{2,1}.d1);
            ig1s = cat(1,xy{1,2}.g,xy{2,1}.g);
        end
        [ir,id,istat]   = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
        coefs(:,cf)     = ir;
        pval(:,cf)      = istat.p;
        r_ci(:,cf,:)    = istat.se*1.96*[-1 1]+ir;

        xy{3,samediff}.y      = id0s;
        xy{3,samediff}.s0     = is0s;
        xy{3,samediff}.s1     = is1s;
        xy{3,samediff}.d1     = kd1s;
        xy{3,samediff}.g      = ig1s;
        xy{3,samediff}.yhat   = glmval(ir,X0(:,2:end),'probit');
        cf = cf + 1;
    end

    npar        = size(r_ci,1);
    pval_fdr    = NaN(npar,6);
    for iseq = 1:6
        [~,~,~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    save('results_GeneralMagnitude.mat','coefs','r_ci','pval_fdr','xy','pval');

    figure(2)

    clf
    colors = lines(4);
    dx = [-1.5 -0.5 0.5 1.5]*0.05;
    for iseq = 1:6
        subplot(1,6,iseq)
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

    %% Motor

    coefs   = [];
    pval    = [];
    r_ci    = [];
    xy      = cell(2,2);
    bis = [1 2];
    for ifeature = [1 2] % {ring, pitch}
        for ih = 1:2
            hs      = [];
            d0s     = [];
            s0s     = [];
            rt0s    = [];
            d01s    = [];
            d02s    = [];
            s1s     = [];
            rt1s    = [];
            g1s     = [];
            jd1s    = [];
            for ig = [1 2] % {2AFC, 4AFC}
                for iSub = 1:nSub(4)
                    par     = pars{4}{iSub};
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

                    y0          = NaN(size(d0));
                    y0(d0>bis(ig))   = 1;
                    y0(d0<=bis(ig))  = 0;

                    [id01,id02] = binarizer(d1,ones(size(d1))*bis(ig)*2);

                    hs      = cat(2,hs,iSameDiffHand);
                    d0s     = cat(2,d0s,y0);
                    s0s     = cat(2,s0s,zscore_HL(s0)); % block-based z-score
                    rt0s    = cat(2,rt0s,zscore_HL(rt0));
                    d01s    = cat(2,d01s,zscore_HL(id01));
                    if ig == 2
                        imat        = NaN(size(id02));
                        iInd        = ~isnan(id02);
                        imat(iInd)  = id02(iInd);
                        imat        = zscore_HL(imat); % block-based z-score
                        id02(iInd)  = imat(iInd);
                    else
                        id02        = zeros(size(id02));
                        id02(isnan(id01)) = NaN;
                    end
                    d02s    = cat(2,d02s,id02);
                    s1s     = cat(2,s1s,zscore_HL(s1));
                    rt1s    = cat(2,rt1s,zscore_HL(rt1));
                    g1s     = cat(2,g1s,ig*ones(size(rt1)));
                    jd1s    = cat(2,jd1s,d1);
                end
            end

            jInd    = hs==ih & ~isnan(d0s + d01s + d02s);

            d0s(~jInd)      = NaN;
            s0s(~jInd)      = NaN;
            rt0s(~jInd)     = NaN;
            d01s(~jInd)     = NaN;
            d02s(~jInd)     = NaN;
            s1s(~jInd)      = NaN;
            rt1s(~jInd)     = NaN;

            ivalid      = ~isnan(sum([d0s(:) s0s(:) rt0s(:) d01s(:) d02s(:) s1s(:) rt1s(:)],2));

            d0s         = d0s(ivalid); % final z-score
            s0s         = zscore(s0s(ivalid));
            rt0s        = zscore(rt0s(ivalid));

            id01s       = zscore(d01s(ivalid));
            uds         = d02s(ivalid);
            iInd        = uds == 0;
            kds         = NaN(size(uds));
            kds(~iInd)  = uds(~iInd);
            kds         = zscore_HL(kds);
            kds(iInd)   = 0;
            id02s       = kds;

            s1s         = zscore(s1s(ivalid));
            rt1s        = zscore(rt1s(ivalid));
            g1s         = g1s(ivalid);
            jd1s        = jd1s(ivalid);

            ids_orig    = id01s;
            ids_cont    = id02s;

            X0_orig         = [d0s ids_orig s0s s1s ids_orig.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) ids_orig.*rt0s s0s.*rt0s s1s.*rt0s ids_orig.*rt1s s1s.*rt1s];
            X0_cont         = [ids_cont ids_cont.*rt0s ids_cont.*rt1s];
            X0              = [X0_orig X0_cont];
            [ir,id,istat]   = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
            coefs(:,2*(ifeature-1)+ih)       = ir;
            pval(:,2*(ifeature-1)+ih)    = istat.p;
            r_ci(:,2*(ifeature-1)+ih,:)  = istat.se*1.96*[-1 1]+ir;

            xy{ifeature,ih}.y      = d0s;
            xy{ifeature,ih}.s0     = s0s;
            xy{ifeature,ih}.s1     = s1s;
            xy{ifeature,ih}.d1     = jd1s;
            xy{ifeature,ih}.g      = g1s;
            xy{ifeature,ih}.yhat   = glmval(ir,X0(:,2:end),'probit');
        end
    end

    npar        = size(r_ci,1);
    h_fdr       = NaN(npar,4);
    pval_fdr    = NaN(npar,4);
    p_thre      = NaN(1,4);
    for iseq = 1:4
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    save('results_Motor.mat','coefs','r_ci','pval_fdr','xy','pval');

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

    %% PreCur

    coefs       = [];
    pval        = [];
    r_ci        = [];
    xy          = cell(1,2);
    bis         = [1 1 1 4];
    bis1        = [1 4 1 1];
    for icomb = 1:2
        d0s     = [];
        s0s     = [];
        rt0s    = [];
        d01s    = [];
        d02s    = [];
        d03s    = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        jd1s    = [];
        for iseq = (1:2)+2*(icomb-1) % 1:2(#)->2(cv), 2:8(#)->2(cv), 3:2(cv)->2(#), 4:2(cv)->8(#)
            for iSub = 1:nSub(5)
                par     = pars{5}{iSub};

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

                y0                  = NaN(size(d0));
                y0(d0>bis(iseq))    = 1;
                y0(d0<=bis(iseq))   = 0;
                y0                  = y0(:);

                [id01,id02,id03] = binarizer(d1,ones(size(d1))*bis1(iseq)*2);

                d0s     = cat(2,d0s,y0);
                s0s     = cat(2,s0s,zscore_HL(s0)); % block-based z-score
                rt0s    = cat(2,rt0s,zscore_HL(rt0));
                d01s    = cat(2,d01s,zscore_HL(id01));
                if iseq ~= 2
                    d02s    = cat(2,d02s,id02);
                    d03s    = cat(2,d03s,id03);
                else
                    d02s    = cat(2,d02s,zscore_HL(id02));
                    d03s    = cat(2,d03s,zscore_HL(id03));
                end
                s1s     = cat(2,s1s,zscore_HL(s1));
                rt1s    = cat(2,rt1s,zscore_HL(rt1));
                g1s     = cat(2,g1s,(2-mod(iseq,2))*ones(size(rt1)));
                jd1s    = cat(2,jd1s,d1);
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) rt0s(:) d01s(:) d02s(:) d03s(:) s1s(:) rt1s(:)],2));

        d0s         = d0s(ivalid);
        s0s         = zscore(s0s(ivalid)); % final z-score
        rt0s        = zscore(rt0s(ivalid));

        uds         = d01s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d1s         = kds;

        uds         = d02s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d2s         = kds;

        uds         = d03s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d3s         = kds;

        d23s        = [d2s d3s];

        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);
        jd1s        = jd1s(ivalid);

        X0_orig         = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt0s s0s.*rt0s s1s.*rt0s d1s.*rt1s s1s.*rt1s];
        X0_cont         = [d23s d23s.*rt0s d23s.*rt1s];

        X0              = [X0_orig X0_cont];
        [ir,id,istat]   = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
        coefs(:,icomb)  = ir;
        pval(:,icomb)    = istat.p;
        r_ci(:,icomb,:)  = istat.se*1.96*[-1 1]+ir;

        xy{icomb}.y      = d0s;
        xy{icomb}.s0     = s0s;
        xy{icomb}.s1     = s1s;
        xy{icomb}.d1     = jd1s;
        xy{icomb}.g      = g1s;
        xy{icomb}.yhat   = glmval(ir,X0(:,2:end),'probit');
    end

    npar        = size(r_ci,1);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pval(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = r_ci(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    save('results_PreCur.mat','coefs','r_ci','pval_fdr','xy','pval');

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

    %% scaling - hard

    coefs   = [];
    pvals   = [];
    CIs     = [];
    xy      = cell(1,2);
    bis     = [1 4];
    for iScaleRepro = 1:2 % 1:scaling, 2:reproduction
        d0s     = [];
        s0s     = [];
        d01s    = [];
        d02s    = [];
        d03s    = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        jd1s    = [];
        for isequence = 1:2 % 1:2->R, 2:8->R, 3:R->2, 4:R->8
            for iSub = 1:nSub(6)
                par     = pars{6}{iSub};

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

                [id01,id02,id03] = binarizer(id1,ones(size(id1))*bis(isequence)*2);

                if iScaleRepro == 1
                    d0s = cat(2,d0s,id0);
                else
                    d0s = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                end
                s0s     = cat(2,s0s,zscore_HL(is0));
                d01s    = cat(2,d01s,zscore_HL(id01));
                if isequence ~= 2
                    d02s    = cat(2,d02s,id02);
                    d03s    = cat(2,d03s,id03);
                else
                    d02s    = cat(2,d02s,zscore_HL(id02));
                    d03s    = cat(2,d03s,zscore_HL(id03));
                end
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));
                jd1s    = cat(2,jd1s,id1);
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d01s(:) d02s(:) d03s(:) s1s(:) rt1s(:)],2));
        if iScaleRepro == 1
            d0s     = d0s(ivalid);
        else
            d0s     = zscore(d0s(ivalid)); % final z-score
        end
        s0s         = zscore(s0s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);
        jd1s        = jd1s(ivalid);

        uds         = d01s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d1s         = kds;

        uds         = d02s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d2s         = kds;

        uds         = d03s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d3s         = kds;

        d23s        = [d2s d3s];
        X0_orig     = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        X0_cont     = [d23s d23s.*rt1s];

        X0          = [X0_orig X0_cont];


        if iScaleRepro == 1
            [ir,id,istat]       = glmfit(X0(:,2:end),X0(:,1),'normal','link','probit'); % normal distribution has much greater fit to binomial distribution
        else
            [ir,~,istat]        = glmfit(X0(:,2:end),X0(:,1));
        end
        coefs(:,iScaleRepro)    = ir;
        pvals(:,iScaleRepro)    = istat.p;
        CIs(:,iScaleRepro,:)    = istat.se*1.96*[-1 1]+ir;

        xy{iScaleRepro}.y      = d0s;
        xy{iScaleRepro}.s0     = s0s;
        xy{iScaleRepro}.s1     = s1s;
        xy{iScaleRepro}.d1     = jd1s;
        xy{iScaleRepro}.g      = g1s;
        if iScaleRepro == 1
            xy{iScaleRepro}.yhat   = glmval(ir,X0(:,2:end),'probit');
        else
            xy{iScaleRepro}.yhat   = glmval(ir,X0(:,2:end),'identity');
        end
    end

    npar        = size(CIs,1);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pvals(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    CIs         = CIs(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    save('results_Scaling.mat','coefs','r_ci','pval_fdr','xy','pvals');

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
            ici     = squeeze(CIs(iInd,iseq,:));
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

    %% scaling - easy

    coefs   = [];
    pvals   = [];
    CIs     = [];
    xy      = cell(1,2);
    for iScaleRepro = 1:2 % 1:scaling, 2:reproduction
        d0s     = [];
        s0s     = [];
        d01s    = [];
        d02s    = [];
        d03s    = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        jd1s    = [];
        for isequence = 1:2 % 1:2->R, 2:8->R, 3:R->2, 4:R->8
            for iSub = 1:nSub(7)
                par     = pars{7}{iSub};

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

                [id01,id02,id03] = binarizer(id1,ones(size(id1))*bis(isequence)*2);

                if iScaleRepro == 1
                    d0s = cat(2,d0s,id0);
                else
                    d0s = cat(2,d0s,zscore_HL(id0)); % block-based z-score
                end
                s0s     = cat(2,s0s,zscore_HL(is0));
                d01s    = cat(2,d01s,zscore_HL(id01));
                if isequence ~= 2
                    d02s    = cat(2,d02s,id02);
                    d03s    = cat(2,d03s,id03);
                else
                    d02s    = cat(2,d02s,zscore_HL(id02));
                    d03s    = cat(2,d03s,zscore_HL(id03));
                end
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                g1s     = cat(2,g1s,isequence*ones(size(irt1)));
                jd1s    = cat(2,jd1s,id1);
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d01s(:) d02s(:) d03s(:) s1s(:) rt1s(:)],2));

        if iScaleRepro == 1
            d0s     = d0s(ivalid);
        else
            d0s     = zscore(d0s(ivalid)); % final z-score
        end
        s0s         = zscore(s0s(ivalid));
        s1s         = zscore(s1s(ivalid));
        rt1s        = zscore(rt1s(ivalid));
        g1s         = g1s(ivalid);
        jd1s        = jd1s(ivalid);

        uds         = d01s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d1s         = kds;

        uds         = d02s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d2s         = kds;

        uds         = d03s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d3s         = kds;

        d23s        = [d2s d3s];
        X0_orig     = [d0s d1s s0s s1s d1s.*(g1s-1) s0s.*(g1s-1) s1s.*(g1s-1) d1s.*rt1s s1s.*rt1s];
        X0_cont     = [d23s d23s.*rt1s];

        X0          = [X0_orig X0_cont];

        if iScaleRepro == 1
            [ir,id,istat]       = glmfit(X0(:,2:end),X0(:,1),'normal','link','probit'); % normal distribution has much greater fit to binomial distribution
            % [ir,id,istat]       = glmfit(X0(:,2:end),X0(:,1),'binomial','link','probit');
        else
            [ir,~,istat]        = glmfit(X0(:,2:end),X0(:,1));
        end
        coefs(:,iScaleRepro)    = ir;
        pvals(:,iScaleRepro)    = istat.p;
        CIs(:,iScaleRepro,:)    = istat.se*1.96*[-1 1]+ir;

        xy{iScaleRepro}.y      = d0s;
        xy{iScaleRepro}.s0     = s0s;
        xy{iScaleRepro}.s1     = s1s;
        xy{iScaleRepro}.d1     = jd1s;
        xy{iScaleRepro}.g      = g1s;
        if iScaleRepro == 1
            xy{iScaleRepro}.yhat   = glmval(ir,X0(:,2:end),'probit');
        else
            xy{iScaleRepro}.yhat   = glmval(ir,X0(:,2:end),'identity');
        end
    end

    npar        = size(CIs,1);
    h_fdr       = NaN(npar,2);
    pval_fdr    = NaN(npar,2);
    p_thre      = NaN(1,2);
    for iseq = 1:2
        [h_fdr(:,iseq),p_thre(iseq),~,pval_fdr(:,iseq)] = fdr_bh(pvals(:,iseq));
    end

    % discard the bias term
    coefs       = coefs(2:end,:);
    r_ci        = CIs(2:end,:,:);
    pval_fdr    = pval_fdr(2:end,:);

    save('results_Scaling_LongStimNoMask.mat','coefs','r_ci','pval_fdr','xy','pvals');

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

for z = 1
    %% ITI varying

    nSub        = nSubs(9);
    inT         = 33;

    pExcs       = 0:0.01:0.3;
    nExc        = length(pExcs);

    results     = [];
    acoefs      = cell(1,3);
    ar_ci       = cell(1,3);
    apval_fdr   = cell(1,3);
    apvals      = [];
    for i = 1:nExc
        pExc        = pExcs(i);
        pminMax     = [pExc 1-pExc];
        clear xy

        d0s     = [];
        s0s     = [];
        rt0s    = [];
        g0s     = [];
        t0s     = [];
        d1s     = [];
        s1s     = [];
        rt1s    = [];
        g1s     = [];
        t1s     = [];
        jd1s    = [];
        d01s    = [];
        d02s    = [];
        d03s    = [];
        for ig = [1 2]
            for iSub = 1:nSub
                iInd    = data{9}.subj==iSub;
                is0     = data{9}.Stm(iInd);
                id0     = data{9}.Chc(iInd);
                irt0    = data{9}.RT(iInd);
                ig0     = (data{9}.granule(iInd)-2)/6+1;
                it0     = (data{9}.iti(iInd)-2)/6+1;
                nR      = length(is0)/inT;

                [~, rank]   = sort(is0, 1, 'ascend'); % 오름차순 정렬
                [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
                ip0         = order/size(is0,1);

                is0     = reshape(is0,[inT,nR]);
                id0     = reshape(id0,[inT,nR]);
                irt0    = reshape(irt0,[inT,nR]);
                ig0     = reshape(ig0,[inT,nR]);
                it0     = reshape(it0,[inT,nR]);
                ip0     = reshape(ip0,[inT,nR]);

                is1     = [NaN(1,nR); is0(1:end-1,:)];
                id1     = [NaN(1,nR); id0(1:end-1,:)];
                irt1    = [NaN(1,nR); irt0(1:end-1,:)];
                ig1     = [NaN(1,nR); ig0(1:end-1,:)];
                it1     = [NaN(1,nR); it0(1:end-1,:)];
                ip1     = [NaN(1,nR); ip0(1:end-1,:)];

                iInd        = (ip0<pminMax(1) | ip0>pminMax(2));

                id0(iInd)   = NaN;
                is0(iInd)   = NaN;
                irt0(iInd)  = NaN;
                ig0(iInd)   = NaN;
                it0(iInd)   = NaN;
                id1(iInd)   = NaN;
                is1(iInd)   = NaN;
                irt1(iInd)  = NaN;
                ig1(iInd)   = NaN;
                it1(iInd)   = NaN;

                iy = NaN(size(id0));
                iy(ig0==1 & id0==1) = 0;
                iy(ig0==1 & id0==2) = 1;
                iy(ig0==2 & id0<=4) = 0;
                iy(ig0==2 & id0>4)  = 1;

                jd1         = NaN(size(id1));
                jd1(ig1==1) = id1((ig1==1));
                [jd01,jd02,jd03] = binarizer(jd1,ones(size(jd1))*2);

                kd1         = NaN(size(id1));
                kd1(ig1==2) = id1((ig1==2));
                [kd01,kd02,kd03] = binarizer(kd1,ones(size(kd1))*8);
                
                jd01        = zscore_HL(jd01);
                kd01        = zscore_HL(kd01);
                kd01(~isnan(jd01)) = jd01(~isnan(jd01));

                kd02        = zscore_HL(kd02);
                kd02(jd02==0) = 0;
                
                kd03        = zscore_HL(kd03);
                kd03(jd03==0) = 0;

                d0s     = cat(2,d0s,iy);
                s0s     = cat(2,s0s,zscore_HL(is0)); % block-based z-score
                rt0s    = cat(2,rt0s,zscore_HL(irt0));
                d1s     = cat(2,d1s,zscore_HL(jd1));
                s1s     = cat(2,s1s,zscore_HL(is1));
                rt1s    = cat(2,rt1s,zscore_HL(irt1));
                jd1s    = cat(2,jd1s,id1);

                d01s    = cat(2,d01s,kd01);
                d02s    = cat(2,d02s,kd02);
                d03s    = cat(2,d03s,kd03);

                g0s     = cat(2,g0s,ig0);
                t0s     = cat(2,t0s,it0);
                g1s     = cat(2,g1s,ig1);
                t1s     = cat(2,t1s,it1);
            end
        end

        ivalid      = ~isnan(sum([d0s(:) s0s(:) d01s(:) d02s(:) d03s(:) rt0s(:) s1s(:) rt1s(:)],2));

        d0s         = d0s(ivalid);
        s0s         = zscore(s0s(ivalid));
        s1s         = zscore(s1s(ivalid));

        uds         = d01s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d1s         = kds;

        uds         = d02s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d2s         = kds;

        uds         = d03s(ivalid);
        iInd        = uds == 0;
        kds         = NaN(size(uds));
        kds(~iInd)  = uds(~iInd);
        kds         = zscore_HL(kds);
        kds(iInd)   = 0;
        d3s         = kds;

        rt0s        = zscore(rt0s(ivalid));
        rt1s        = zscore(rt1s(ivalid));

        g0s         = g0s(ivalid);
        g1s         = g1s(ivalid);

        t0s         = t0s(ivalid);
        t1s         = t1s(ivalid);

        jd1s        = jd1s(ivalid);


        ix          = [d1s s0s s1s];
        d23s        = [d2s d3s];
        X0_orig     = [d0s ix ix.*(g1s-1) ix.*(t1s-3/2) ix.*(g1s-1).*(t1s-3/2) ix.*rt0s [d1s s1s].*rt1s];
        X0_cont     = [d23s d23s.*(t1s-3/2) d23s.*rt0s d23s.*rt1s];
        X00         = [X0_orig X0_cont];
        [ir,id,istat]   = glmfit(X00(:,2:end),X00(:,1),'binomial','link','probit');

        X0_orig                 = [d0s ix ix.*(g1s-1) ix.*rt0s [d1s s1s].*rt1s];
        X0_cont                 = [d23s d23s.*rt0s d23s.*rt1s];
        X0                      = [X0_orig X0_cont];
        [ir2_1,~,istat2_1]      = glmfit(zscore(X0(t1s==1,2:end)),X0(t1s==1,1),'binomial','link','probit');
        [ir2_2,~,istat2_2]      = glmfit(zscore(X0(t1s==2,2:end)),X0(t1s==2,1),'binomial','link','probit');

        X0_orig        = [d0s ix ix.*rt0s [d1s s1s].*rt1s];
        X0_cont        = [d23s d23s.*rt0s d23s.*rt1s];
        X0             = [X0_orig X0_cont];
        [~,~,istat1]   = glmfit(zscore(X0(t1s==1 & g1s==1,2:end)),X0(t1s==1 & g1s==1,1),'binomial','link','probit');
        [~,~,istat2]   = glmfit(zscore(X0(t1s==1 & g1s==2,2:end)),X0(t1s==1 & g1s==2,1),'binomial','link','probit');
        [~,~,istat3]   = glmfit(zscore(X0(t1s==2 & g1s==1,2:end)),X0(t1s==2 & g1s==1,1),'binomial','link','probit');
        [~,~,istat4]   = glmfit(zscore(X0(t1s==2 & g1s==2,2:end)),X0(t1s==2 & g1s==2,1),'binomial','link','probit');

        coefs   = ir;
        pvals   = istat.p;
        r_ci    = istat.se*1.96*[-1 1]+ir;

        xy.y      = d0s;
        xy.s0     = s0s;
        xy.s1     = s1s;
        xy.d1     = jd1s;
        xy.g      = g1s;
        xy.t0     = t0s;
        xy.t1     = t1s;
        xy.yhat   = glmval(ir,X00(:,2:end),'probit');

        [~,~,~,pval_fdr] = fdr_bh(pvals);
        [~,~,~,pval_fdr2_1] = fdr_bh(istat2_1.p);
        [~,~,~,pval_fdr2_2] = fdr_bh(istat2_2.p);

        % discard the bias term
        coefs       = coefs(2:end);
        r_ci        = r_ci(2:end,:);
        pval_fdr    = pval_fdr(2:end);

        ir2_1       = ir2_1(2:end);
        r_ci2_1     = istat2_1.se(2:end)*1.96*[-1 1]+ir2_1;
        pval_fdr2_1 = pval_fdr2_1(2:end);

        ir2_2       = ir2_2(2:end);
        r_ci2_2     = istat2_2.se(2:end)*1.96*[-1 1]+ir2_2;
        pval_fdr2_2 = pval_fdr2_2(2:end);

        results(i).pExc = pExc;
        results(i).coefs = coefs;
        results(i).r_ci = r_ci;
        results(i).pval_fdr = pval_fdr;
        results(i).xy = xy;
        results(i).pvals = pvals;
        results(i).stats{1} = istat1;
        results(i).stats{2} = istat2;
        results(i).stats{3} = istat3;
        results(i).stats{4} = istat4;
        results(i).stats{5} = istat2_1;
        results(i).stats{6} = istat2_2;

        acoefs{1}      = cat(2,acoefs{1},coefs);
        ar_ci{1}       = cat(3,ar_ci{1},r_ci);
        apval_fdr{1}   = cat(2,apval_fdr{1},pval_fdr);

        acoefs{2}      = cat(2,acoefs{2},ir2_1);
        ar_ci{2}       = cat(3,ar_ci{2},r_ci2_1);
        apval_fdr{2}   = cat(2,apval_fdr{2},pval_fdr2_1);

        acoefs{3}      = cat(2,acoefs{3},ir2_2);
        ar_ci{3}       = cat(3,ar_ci{3},r_ci2_2);
        apval_fdr{3}   = cat(2,apval_fdr{3},pval_fdr2_2);

        apvals      = cat(2,apvals,pvals(2:end));
    end
    results(end+1).pExcs = pExcs;

    save('results_ITI.mat','results','acoefs','ar_ci','apval_fdr','apvals');
end