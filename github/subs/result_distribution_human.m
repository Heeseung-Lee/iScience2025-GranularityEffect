clear all; clc


dthetabin       = 0.3;
drho            = 0.1;

tbin            = (-3+dthetabin/2):dthetabin:(3-dthetabin/2);
rangerho        = (0+drho/2):drho:(1-drho/2);

nbin            = 8;
pbin            = linspace(0,1,nbin+1);

[pars, nSub]    = RingGaborDataLoad;

ifeature    = 2;
cg          = 1;
Gs          = [2 4 8];
dists       = cell(1,2);
coefs       = cell(1,2);
for ig = [1 2 3]

    inG         = Gs(ig);
    dists{cg}   = NaN(nbin,inG,nbin,inG,nSub); % [prechoice, prestim, curchoice, subj]
    for iSub = 1:nSub
        par     = pars{iSub,ifeature};
        iIndR   = find((par.condition == ig) & (par.StairTrainTest==3));
        d0      = par.Chc(:,iIndR);
        s0      = par.Stm(:,iIndR);
        inT     = size(s0,1);
        inR     = size(s0,2);
        
        xbin        = icdf('normal',pbin,par.iref,par.isig);
        i0          = NaN(inT,inR);
        for ibin = 1:nbin
            iInd        = s0>=xbin(ibin) & s0<xbin(1+ibin);
            i0(iInd)    = ibin;
        end

        i1  = [NaN(1,size(i0,2)); i0(1:end-1,:)];
        d1  = [NaN(1,size(i0,2)); d0(1:end-1,:)];
        s1  = [NaN(1,size(i0,2)); s0(1:end-1,:)];

        for is1 = 1:nbin
            for id1 = 1:inG
                for id0 = 1:inG
                    for is0 = 1:nbin
                        iInd                                = i1==is1 & d1 == id1 & d0 == id0 & i0 == is0;
                        dists{cg}(is0,id1,is1,id0,iSub)     = sum(iInd(:));
                    end
                end
            end
        end

        X                   = [d0(:) s0(:) d1(:) s1(:)];
        X                   = X(~isnan(sum(X,2)),:);
        
        % multinomial logistic regression
        ix      = table(discretize(-X(:,1),[(-inG:-1)-0.5 -0.5],'categorical'),VariableNames="d0");
        X       = zscore(X(:,2:end));
        X       = cell2table(num2cell(X),'VariableNames',{'i0','d1','i1'});
        X       = [ix X];
        imodel  = fitmnr(X,'d0',ModelType="ordinal");
        icoef   = imodel.Coefficients.Value;

        % logistic regression
        % X(:,1)              = X(:,1)>Gs(ig)/2;
        % X(:,2:end)          = zscore(X(:,2:end));
        % icoef               = glmfit(X(:,2:end),X(:,1),'binomial','link','probit');
       
        % linear regression
        % X                   = zscore(X);
        % icoef               = glmfit(X(:,2:end),X(:,1));

        coefs{ig}(iSub,:)   = icoef(end-2:end);
    end
    cg = cg + 1;
end

% effect of s0

nGs = [2 4 8];

figure(1+4*(ifeature-1))
clf
colormap(copper)
colors = lines(3);
for ig = 1:3
    inG = nGs(ig);
    pdfs = [];
    for id1 = 1:inG
        for is1 = 1:nbin
            idist   = squeeze(mean(dists{ig}(:,id1,is1,:,:),5,'omitnan')); %[d1, d0]
            ipdf    = idist./sum(idist,2);
            pdfs    = cat(3,pdfs,ipdf);
        end
    end

    subplot(3,2,2*(ig-1)+1);
    hold on
    ipdf    = mean(pdfs,3,'omitnan')';
    h       = imagesc(ipdf);
    set(h, 'AlphaData', 1-isnan(ipdf))
    xlim([0.5 0.5+size(ipdf,2)])
    ylim([0.5 0.5+size(ipdf,1)])
    set(gca,'ytick',1:inG)
    xlabel('current stimulus')
    ylabel('current choice')
    clim([0 1])


    subplot(3,2,2*ig);
    hold on
    ipdf    = [sum(ipdf(1:inG/2,:),1,'omitnan'); sum(ipdf((inG/2+1):end,:),1,'omitnan')];
    dpdf = ipdf(2,:);
    plot(dpdf,'-o','linewidth',2,'color',colors(ig,:))
    plot([1 nbin],[0.5 0.5],'k--')
    ylabel('current choice = large')
    xlabel('current stimulus')
    ylim([0 1])
    xlim([1 nbin])
end

%
%
% effect of d1

figure(2+4*(ifeature-1))
clf
colors = lines(3);
colormap(copper)
for ig = 1:3
    pdfs = [];
    for is1 = 1:nbin
        for is0 = 1:nbin
            idist   = squeeze(mean(dists{ig}(is0,:,is1,:,:),5,'omitnan')); %[d1, d0]
            ipdf    = idist./sum(idist,2);
            pdfs    = cat(3,pdfs,ipdf);
        end
    end
    inG     = nGs(ig);

    subplot(3,2,2*(ig-1)+1);
    hold on
    ipdf    = mean(pdfs,3,'omitnan')';
    h       = imagesc(ipdf);
    set(h, 'AlphaData', 1-isnan(ipdf))
    xlim([0.5 0.5+size(ipdf,1)])
    ylim([0.5 0.5+size(ipdf,1)])
    set(gca,'xtick',1:inG,'ytick',1:inG)
    xlabel('previous choice')
    ylabel('current choice')
    clim([0 1])

    subplot(3,2,2*ig);
    hold on
    ipdf    = [sum(ipdf(1:inG/2,:),1,'omitnan'); sum(ipdf((inG/2+1):end,:),1,'omitnan')];
    dpdf = ipdf(2,:);
    plot(dpdf,'-o','linewidth',2,'color',colors(ig,:))
    ylabel('current choice = large')
    xlabel('previous choice')
    plot([1 inG],[0.5 0.5],'k--')
    ylim([0 1])
    xlim([1 inG])
end

% effect of s1

figure(3+4*(ifeature-1))
clf
colormap(copper)
colors = lines(3);
for ig = 1:3
    inG = nGs(ig);
    pdfs = [];
    for id1 = 1:inG
        for is0 = 1:nbin
            idist   = squeeze(mean(dists{ig}(is0,id1,:,:,:),5,'omitnan')); %[d1, d0]
            ipdf    = idist./sum(idist,2);
            pdfs    = cat(3,pdfs,ipdf);
        end
    end

    subplot(3,2,2*(ig-1)+1);
    hold on
    ipdf    = mean(pdfs,3,'omitnan')';
    h       = imagesc(ipdf);
    set(h, 'AlphaData', 1-isnan(ipdf))
    xlim([0.5 0.5+size(ipdf,2)])
    ylim([0.5 0.5+size(ipdf,1)])
    set(gca,'ytick',1:inG)
    xlabel('previous stimulus')
    ylabel('current choice')
    clim([0 1])

    subplot(3,2,2*ig);
    hold on
    ipdf    = [sum(ipdf(1:inG/2,:),1,'omitnan'); sum(ipdf((inG/2+1):end,:),1,'omitnan')];
    dpdf = ipdf(2,:);
    plot(dpdf,'-o','linewidth',2,'color',colors(ig,:))
    ylabel('current choice = large')
    xlabel('previous stimulus')
    plot([1 nbin],[0.5 0.5],'k--')
    ylim([0 1])
    xlim([1 nbin])

end

%

ncoef = 3;
figure(4*ifeature)
clf
dx= 0.2;
names = {'s0','c1','s1'};
for icoef = 1:ncoef
    subplot(1,3,icoef)
    hold on
    jcoef       = [coefs{1}(:,icoef) coefs{2}(:,icoef) coefs{3}(:,icoef)];
    m           = mean(jcoef);
    [~,~,ci]    = ttest(jcoef);
    plot(1:3,jcoef,'color',[1 1 1]*0.7)
    for ig = 1:3
        bar(ig,m(ig),'EdgeColor',colors(ig,:),'facecolor','none','linewidth',2)
        plot([ig ig],ci(:,ig),'_-','color',colors(ig,:),'linewidth',2)
    end
    dy = mean(ci(2,:) - ci(1,:))/2;
    p = NaN(3,3);
    for ig = 1:2
        jg          = ig+1;
        [~,ip]      = ttest(jcoef(:,ig),jcoef(:,jg));
        p(ig,jg)    = ip;
        if ip < 0.05
            if ip < 0.001
                ix = [-1 0 1]*dx;
            elseif ip < 0.01
                ix = [-1 1]*dx/2;
            else
                ix = 0;
            end
            plot(ig+0.5+ix,max(ci(:))+1.5*dy,'k','marker','*','linewidth',1.5,'markersize',10)
        end
    end
    icoefs = NaN(nSub,1);
    for iSub = 1:nSub
        kcoef   = jcoef(iSub,:);
        kcoef   = glmfit(1:3,zscore(kcoef),'normal','Link','identity');
        icoefs(iSub) = kcoef(2);
    end
    [~,p] = ttest(icoefs);    
    title(sprintf([names{icoef} ' | p=%.3f'],p))
    set(gca,'xtick',1:3,'ytick',-1.5:0.5:3,'xticklabel',[2 4 8])
    xlabel('G')
    ylabel('coefficient')
end

%%

