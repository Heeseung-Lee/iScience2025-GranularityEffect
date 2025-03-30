% Model 1: only theta-update                {sig_m, dif_t}
% Model 2: only rho-update                  {sig_m, dif_r}
% Model 3: theta & rho-update               {sig_m, dif_t, dif_r}
% Model 4: theta-update & rho-condition     {sig_m, dif_t, dif_r}

% Cond 1: G0=scaling & G1=2
% Cond 2: G0=scaling & G1=8
% Cond 3: G0=reproduction & G1=2
% Cond 4: G0=reproduction & G1=8
% Cond 5: G0=2 & G1=2
% Cond 6: G0=2 & G1=8
% Cond 7: G0=8 & G1=2

clear all; clc
addpath('/Volumes/CSNL_new/people/HSL/projects/granularity/codes/subs')
dir0        = '/Volumes/CSNL_new/people/HSL/projects/granularity/results/Figures/models';

nIter_total = 100;
nIter_plot  = 16;
nIter_data  = 1000;

neach       = 10;
sig_m       = linspace(0.4,0.8,neach); % linspace(0.4,0.7,4)
dif_t       = linspace(0.3,0.7,neach); % linspace(0.2,0.5,4)
dif_r       = linspace(0.2,0.6,neach); % linspace(0.2,0.4,4)

models      = 1:4;
nmodel      = length(models);
nT          = 2000;
N_renew     = 5;
nm          = length(sig_m);
nt          = length(dif_t);
nr          = length(dif_r);

pCI_plot            = 95;
Gs                  = [2 8];
nG                  = length(Gs);
rhobounds           = cell(1,nG);
for ig = 1:nG
    rhobounds{ig}   = linspace(0,1,Gs(ig)+1);
end

dtheta              = 0.01;
drho                = 0.01;
rangetheta          = (-7+dtheta/2):dtheta:(7-dtheta/2);
rangerho            = (0+drho/2):drho:(1-drho/2);
rangerho0           = (-4+drho/2):drho:(5-drho/2);
nrho                = length(rangerho);

%
savedir     = sprintf([dir0 '/sig-m%.1f%.1f dif-t%.1f%.1f dif-r%.1f%.1f nEach%d nT%d renew%d'],sig_m(1),sig_m(end),dif_t(1),dif_t(end),dif_r(1),dif_r(end),neach,nT,N_renew);
if isempty(dir(savedir))
    mkdir(savedir)
end
cd(savedir)

for iIter = 1:nIter_total
    if isempty(dir([savedir '/Iter' num2str(iIter) '.mat'])) && isempty(dir([savedir '/Iter' num2str(iIter) '_compute.mat'])) && iIter <= nIter_plot
        save([savedir '/Iter' num2str(iIter) '_compute.mat'],'iIter')

        betas = cell(1,max(models));
        param = cell(1,max(models));
        for imodel = models
            switch imodel
                case 1 % {sig_m, dif_t}
                    npar = nm*nt;
                    pars = NaN(npar,2);
                    c = 1;
                    for im = 1:nm
                        for it = 1:nt
                            pars(c,:)   = [sig_m(im) dif_t(it)];
                            c = c + 1;
                        end
                    end
                case 2 % {sig_m, dif_r}
                    npar = nm*nr;
                    pars = NaN(npar,2);
                    c = 1;
                    for im = 1:nm
                        for ir = 1:nr
                            pars(c,:)   = [sig_m(im) dif_r(ir)];
                            c = c + 1;
                        end
                    end
                case {3,4} % {sig_m, dif_t, dif_r}
                    npar = nm*nt*nr;
                    pars = NaN(npar,3);
                    c = 1;
                    for im = 1:nm
                        for ir = 1:nr
                            for it = 1:nt
                                pars(c,:)   = [sig_m(im) dif_t(it) dif_r(ir)];
                                c = c + 1;
                            end
                        end
                    end
            end
            param{imodel} = pars;
            betas{imodel} = NaN(3,7,npar);

            for ipar = 1:npar
                fprintf('iIter=%d/%d | imodel=%d/%d | ipar=%d/%d\n',iIter,nIter_total,imodel,nmodel,ipar,npar)

                isig_m = pars(ipar,1);
                switch imodel
                    case 1
                        idif_t = pars(ipar,2);
                    case 2
                        idif_r = pars(ipar,2);
                    case {3,4}
                        idif_t = pars(ipar,2);
                        idif_r = pars(ipar,3);

                        if imodel == 4
                            prior_rho = cell(1,max(Gs));
                            for iG = Gs
                                iprior_rho  = NaN(iG,nrho);
                                for iclass = 1:iG
                                    iprior_rho(iclass,:)    = ConditionalPrior_Rho(iG,iclass,idif_r,rangerho0);
                                end
                                prior_rho{iG}   = iprior_rho;
                            end
                        end
                end

                theta       = normrnd(0,1,nT,1);
                m           = normrnd(theta,isig_m,nT,1);

                indG        = randi(2,[nT,1]);
                indScaling  = NaN(nT,1);
                theta_hat   = NaN(nT,1);
                classes     = NaN(nT,1);
                rho_hat     = NaN(nT,1);
                for iT = 1:nT
                    indScaling(iT)  = mod(iT,N_renew);
                    if indScaling(iT) == 1 % initialize prior for every Nth trial
                        imu_pri     = 0;
                        isig_pri    = 1; % uniform
                        irho_pri    = ones(1,nrho)/drho/nrho;
                    end
                    im              = m(iT);

                    % theta inference
                    imu_pos         = (im*isig_pri^2 + imu_pri*isig_m^2)/(isig_pri^2 + isig_m^2);
                    isig_pos        = isig_pri*isig_m/sqrt(isig_pri^2 + isig_m^2);
                    theta_hat(iT)   = imu_pos;

                    % rho inference
                    rhox            = normcdf(rangetheta,imu_pri,isig_pri);
                    jrho_lik        = normpdf(rangetheta,im,isig_m);
                    [~,iInd]        = unique(rhox);
                    irho_lik        = interp1(rhox(iInd),jrho_lik(iInd),rangerho);
                    irho_lik(isnan(irho_lik)) = 0;
                    irho_lik        = irho_lik/sum(irho_lik*drho);
                    irho_pos        = irho_lik.*irho_pri;
                    irho_pos        = irho_pos/sum(irho_pos*drho);
                    irho_hat        = sum(rangerho.*irho_pos*drho);
                    rho_hat(iT)     = irho_hat;
                    %
                    iG  = indG(iT);
                    inG = Gs(iG);
                    for ig = 1:inG
                        irange              = rhobounds{iG}(ig:ig+1);
                        if irho_hat >= irange(1) && irho_hat < irange(2)
                            classes(iT)     = ig;
                            break
                        end
                    end

                    if indScaling(iT) ~= 0 % does not update prior on the testing trial
                        % Module: theta-prior update
                        if imodel == 1 || imodel == 3 || imodel == 4
                            imu_pri   = imu_pos;
                            isig_pri  = sqrt(isig_pos^2 + idif_t^2);
                        end

                        % Module: rho prior-update
                        if imodel == 2 || imodel == 3
                            ix0         = round(rangerho0,5);
                            ix          = round(rangerho,5);
                            iInd        = find(ix0==ix(1)):find(ix0==ix(end));
                            iy          = zeros(length(ix0),1);
                            iy(iInd)    = irho_pos;
                            dx          = ix0(2) - ix0(1);
                            icumy       = cumsum(iy*dx);
                            iInd        = find(icumy>0.0001,1,'first'):find(icumy<0.9999,1,'last');
                            if isempty(iInd)
                                iInd    = find(icumy>0.9999,1,'first');
                            end
                            nx          = length(ix0);
                            ip          = zeros(1,nx);
                            nInd        = length(iInd);
                            for i = 1:nInd
                                kInd    = iInd(i);
                                jx      = ix0 - ix0(kInd);
                                ip      = ip + iy(kInd)*normpdf(jx,0,idif_r);
                            end
                            ip          = ip/sum(ip*dx);
                            ip          = range_truncate_from0to1(ip,ix0);
                            irho_pri    = ip;
                        end

                        % Module: rho prior-condition
                        if imodel == 4
                            irho_pri    = prior_rho{inG}(classes(iT),:);
                        end
                    end

                end

                %
                indG1   = [NaN(1,1); indG(1:end-1)];
                for ig = 1:7
                    switch ig
                        case {1,2,3,4}
                            iInd    = find(indScaling==0 & indG1==ig - 2*(ig>2) - 2*(ig>4));
                        case 5
                            iInd    = find(indScaling==0 & indG==1 & indG1==1);
                        case 6
                            iInd    = find(indScaling==0 & indG==1 & indG1==2);
                        case 7
                            iInd    = find(indScaling==0 & indG==2 & indG1==1);
                    end
                    iInd    = iInd(2:end);
                    if ig < 3
                        X   = [rho_hat(iInd) classes(iInd-1) theta(iInd) theta(iInd-1)];
                    elseif ig < 5
                        X   = [theta_hat(iInd) classes(iInd-1) theta(iInd) theta(iInd-1)];
                    else
                        X   = [classes(iInd) classes(iInd-1) theta(iInd) theta(iInd-1)];
                    end
                    iInd    = ~isnan(sum(X,2));
                    X       = zscore(X(iInd,:));
                    ibeta   = glmfit(X(:,2:end),X(:,1));
                    betas{imodel}(:,ig,ipar) = ibeta(2:end);
                end
            end
        end
        save([savedir '/Iter' num2str(iIter) '.mat'],'betas','param')
        delete([savedir '/Iter' num2str(iIter) '_compute.mat'])
    end
end

%%

completes = NaN(nIter_plot,1);
for iIter = 1:nIter_plot
    completes(iIter,1) = ~isempty(dir([savedir '/Iter' num2str(iIter) '.mat'])) & isempty(dir([savedir '/Iter' num2str(iIter) '_compute.mat']));
end
condname = {'scale','reprod','classG0=2','classG1=2'};

if sum(completes) == nIter_plot

    % simulation load
    betas = cell(1,max(models));
    param = cell(1,max(models));
    for iIter = 1:nIter_plot
        for imodel = models
            iresult = load([savedir '/Iter' num2str(iIter) '.mat']);
            betas{imodel}       = cat(4,betas{imodel},iresult.betas{imodel});
            if iIter == 1
                param{imodel}   = cat(4,param{imodel},iresult.param{imodel});
            end
        end
    end

    %% observation load
    CIs         = cell(1,4);
    load_dir    = ['/Volumes/CSNL_new/people/HSL/projects/granularity/results/Figures/humans/nIter' num2str(nIter_data)];
    load([load_dir '/results_Scaling.mat'],'coefs_bst','coefs');
    for iseq = 1:2
        for ig = 1:2
            icoefs  = squeeze(coefs_bst(:,ig+2*(iseq-1),:));
            ici     = prctile(icoefs,[0 100]+[1 -1]*(100-pCI_plot)/2,2);
            CIs{iseq}(:,:,ig)   = ici;
        end
    end
    load([load_dir '/results_PreCur.mat'],'coefs_bst');
    for iseq = 1:2
        for ig = 1:2
            icoefs  = squeeze(coefs_bst(:,ig+2*(iseq-1),:));
            ici     = prctile(icoefs,[0 100]+[1 -1]*(100-pCI_plot)/2,2);
            CIs{2+iseq}(:,:,ig)     = ici;
        end
    end

    for imodel = models

        mbeta       = mean(betas{imodel},4);
        ipar        = param{imodel};
        npar        = size(ipar,2);
        ncomb       = nchoosek(npar,2);
        coefname    = {'d1','s0','s1'};
        if imodel ~= 2
            parname     = {'sig-m','dif-t','dif-r'};
        else
            parname     = {'sig-m','dif-r'};
        end

        % constraint: data performance at G0=2 classification
        iCI         = mean(CIs{3}(2,:,:),3); 
        dbeta       = squeeze(mbeta(2,5,:));
        IndValid    = cell(ncomb,1);
        cp = 1;
        for ip = 1:npar
            for jp = (ip+1):npar
                ips     = unique(ipar(:,ip));
                jps     = unique(ipar(:,jp));
                nip     = length(ips);
                njp     = length(jps);
                mdbeta  = NaN(nip,njp);
                for i = 1:nip
                    for j = 1:njp
                        iInd        = ipar(:,ip) == ips(i) & ipar(:,jp) == jps(j);
                        mdbeta(i,j) = mean(dbeta(iInd));
                    end
                end
                iValid          = find(iCI(1) < mdbeta & iCI(2) > mdbeta);
                [iInd,jInd]     = ind2sub(size(mdbeta),iValid);
                IndValid{cp}    = [iInd jInd];
                cp = cp + 1;
            end
        end

        for icond = 1:4

            switch icond
                case 1 % scaling G1=2,8
                    dbeta = squeeze(mbeta(:,1,:))';
                case 2 % reproduction G1=2,8
                    dbeta = squeeze(mbeta(:,3,:))';
                case 3 % classification G1=2,8 & G0=2
                    dbeta = squeeze(mbeta(:,5,:))';
                case 4 % classification G1=2 & G0=2,8
                    dbeta = squeeze(mbeta(:,5,:))';
            end

            nbeta   = size(dbeta,2);

            figure(icond+4)
            clf
            for ib = [2 1 3]
                cp = 1;
                for ip = 1:npar
                    for jp = (ip+1):npar
                        subplot(nbeta,ncomb,cp + ncomb*(ib-1))
                        hold on
                        ips = unique(ipar(:,ip));
                        jps = unique(ipar(:,jp));

                        nip = length(ips);
                        njp = length(jps);

                        mdbeta = NaN(nip,njp);
                        for i = 1:nip
                            for j = 1:njp
                                iInd        = ipar(:,ip) == ips(i) & ipar(:,jp) == jps(j);
                                mdbeta(i,j) = mean(dbeta(iInd,ib));
                            end
                        end
                        imagesc(mdbeta)
                        
                        iIndValid   = IndValid{cp};
                        nValid      = size(iIndValid,1);
                        for ival = 1:nValid
                            plot(iIndValid(ival,1),iIndValid(ival,2),'k*')
                        end
                        colorbar
                        title(coefname{ib})
                        xlim([0.5 njp+0.5])
                        ylim([0.5 nip+0.5])
                        if ib == 2
                            clim([0 1])
                        else
                            clim([-0.7 0.7])
                        end
                        cp = cp + 1;
                        xlabel(parname{jp})
                        ylabel(parname{ip})
                        set(gca,'xtick',1:njp,'xticklabel',round(jps,1),'ytick',1:nip,'yticklabel',round(ips,1))
                    end
                end
            end
            saveas(gcf,[num2str(icond) condname{icond} '_G2_M' num2str(imodel) '.png'])


            switch icond
                case 1 % scaling G1=2,8
                    dbeta = squeeze(mbeta(:,2,:) - mbeta(:,1,:))';
                case 2 % reproduction G1=2,8
                    dbeta = squeeze(mbeta(:,4,:) - mbeta(:,3,:))';
                case 3 % classification G1=2,8 & G0=2
                    dbeta = squeeze(mbeta(:,6,:) - mbeta(:,5,:))';
                case 4 % classification G1=2 & G0=2,8
                    dbeta = squeeze(mbeta(:,7,:) - mbeta(:,5,:))';
            end

            figure(icond)
            clf
            for ib = 1:nbeta
                cp = 1;
                for ip = 1:npar
                    for jp = (ip+1):npar
                        subplot(nbeta,ncomb,cp + ncomb*(ib-1))
                        hold on
                        ips = unique(ipar(:,ip));
                        jps = unique(ipar(:,jp));

                        nip = length(ips);
                        njp = length(jps);

                        mdbeta = NaN(nip,njp);
                        for i = 1:nip
                            for j = 1:njp
                                iInd        = ipar(:,ip) == ips(i) & ipar(:,jp) == jps(j);
                                mdbeta(i,j) = mean(dbeta(iInd,ib));
                            end
                        end
                        imagesc(mdbeta)
                        iIndValid   = IndValid{cp};
                        nValid      = size(iIndValid,1);
                        for ival = 1:nValid
                            plot(iIndValid(ival,1),iIndValid(ival,2),'k*')
                        end
                        colorbar
                        title(coefname{ib})
                        xlim([0.5 njp+0.5])
                        ylim([0.5 nip+0.5])
                        clim([-0.2 0.2])
                        cp = cp + 1;
                        xlabel(parname{jp})
                        ylabel(parname{ip})
                        set(gca,'xtick',1:njp,'xticklabel',round(jps,1),'ytick',1:nip,'yticklabel',round(ips,1))
                    end
                end
            end

            saveas(gcf,[num2str(icond) condname{icond} '_Del28_M' num2str(imodel) '.png'])
        end
    end

end





