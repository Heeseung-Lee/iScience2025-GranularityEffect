clear all; clc

dir0 = '/Users/heeseunglee/Documents/project/granularity/github';
addpath([dir0 '/subs'])

nIter       = 100;
nT          = 1000;
npar        = 10;
dr          = 0.01;
rrange0     = (-4+dr/2):dr:(5-dr/2);

for iIter = 1:nIter

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% General settings %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rrange      = (dr/2):dr:(1-dr/2);
    nrr         = length(rrange);
    Gs          = [2 8];
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

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Standard settings %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig_ms      = linspace(0.1,1,npar);
    sig_mms     = linspace(0.4,4,npar);
    sig_pri0s   = linspace(0.4,4,npar);
    sig_rs      = linspace(0.02,0.2,npar);
    dif_rs      = linspace(0.02,0.2,npar);

    pars        = allcomb(sig_ms,sig_mms,sig_pri0s,sig_rs,dif_rs);
    npars       = size(pars,1);

    isavename   = sprintf([dir0 '/result/models/Standard/'...
        'ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_srs%.2f%.2f_drs%.2f%.2f_nT%d_npar%d'],...
        sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),sig_rs(1),sig_rs(end),dif_rs(1),dif_rs(end),nT,npar);
    if isempty(dir(isavename))
        mkdir(isavename)
    end
    ifilename = [isavename '/Iter' num2str(iIter) '.mat'];

    if isempty(dir(ifilename))

        coefs = NaN(8,3,npars); % ncoef, nvar, npar
        perfs = NaN(npars,1);
        for ipar = 1:npars
            [t, indG, thetas, rhos, classes, iperfs] = model_Standard(pars(ipar,:), Gs, nT, nrr, rrange0, rrange, rhobounds, mu_pri0);
            coefs(:,:,ipar) = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
            perfs(ipar) = iperfs;
            fprintf('iIter=%3d/%3d | Standard | ipar=%5d/%5d | iperf=%4.1f\n',iIter,nIter,ipar,npars,perfs(ipar))
        end

        dataStruct = struct("coefs",coefs,"pars",pars,"perfs",perfs);
        save(ifilename,'-fromstruct',dataStruct)

    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Adaptation settings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig_ms      = linspace(0.1,1,npar);
    gains       = linspace(0.1,1,npar);
    sig_pri0s   = linspace(0.4,4,npar);

    pars        = allcomb(sig_ms,gains,sig_pri0s);
    npars       = size(pars,1);

    isavename   = sprintf([dir0 '/result/models/Adaptation/ms%.2f%.2f_gains%.2f%.2f_pri0%.2f%.2f_nT%d_npar%d'],...
        sig_ms(1),sig_ms(end),gains(1),gains(end),sig_pri0s(1),sig_pri0s(end),nT,npar);
    if isempty(dir(isavename))
        mkdir(isavename)
    end
    ifilename = [isavename '/Iter' num2str(iIter) '.mat'];

    if isempty(dir(ifilename))

        coefs = NaN(8,3,npars); % ncoef, nbin, nvar, npar
        perfs = NaN(npars,1);
        for ipar = 1:npars
            [t, indG, thetas, rhos, classes, iperf] = model_Adaptation(pars(ipar,:), Gs, nT, trange, dt, rhobounds, mu_pri0);
            coefs(:,:,ipar) = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
            perfs(ipar) = iperf;

            fprintf('iIter=%3d/%3d | Adaptation | ipar=%5d/%5d | iperf=%4.1f\n',iIter,nIter,ipar,npars,perfs(ipar))
        end

        dataStruct = struct("coefs",coefs,"pars",pars,"perfs",perfs);
        save(ifilename,'-fromstruct',dataStruct)

    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ThetaUpdate settings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig_ms      = linspace(0.1,1,npar);
    sig_mms     = linspace(0.4,4,npar);
    sig_pri0s   = linspace(0.4,4,npar);

    pars        = allcomb(sig_ms,sig_mms,sig_pri0s);
    npars       = size(pars,1);

    isavename   = sprintf([dir0 '/result/models/ThetaUpdate/ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_nT%d_npar%d'],...
        sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),nT,npar);
    if isempty(dir(isavename))
        mkdir(isavename)
    end
    ifilename = [isavename '/Iter' num2str(iIter) '.mat'];

    if isempty(dir(ifilename))

        coefs = NaN(8,3,npars); % ncoef, nbin, nvar, npar
        perfs = NaN(npars,1);
        for ipar = 1:npars
            [t, indG, thetas, rhos, classes, iperf] = model_ThetaUpdate(pars(ipar,:), Gs, nT, rhobounds, mu_pri0);
            coefs(:,:,ipar) = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
            perfs(ipar) = iperf;

            fprintf('iIter=%3d/%3d | Theta Update | ipar=%5d/%5d | iperf=%4.1f\n',iIter,nIter,ipar,npars,perfs(ipar))
        end

        dataStruct = struct("coefs",coefs,"pars",pars,"perfs",perfs);
        save(ifilename,'-fromstruct',dataStruct)

    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ThetaConded settings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig_ms      = linspace(0.1,1,npar);
    sig_mms     = linspace(0.4,4,npar);
    sig_pri0s   = linspace(0.4,4,npar);
    dif_ts      = linspace(0.2,2,npar);

    pars        = allcomb(sig_ms,sig_mms,sig_pri0s,dif_ts);
    npars       = size(pars,1);

    isavename   = sprintf([dir0 '/result/models/ThetaConditioned/ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_ts%.2f%.2f_nT%d_npar%d'],...
        sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),dif_ts(1),dif_ts(end),nT,npar);
    if isempty(dir(isavename))
        mkdir(isavename)
    end
    ifilename = [isavename '/Iter' num2str(iIter) '.mat'];

    if isempty(dir(ifilename))

        coefs = NaN(8,3,npars); % ncoef, nbin, nvar, npar
        perfs = NaN(npars,1);
        for ipar = 1:npars
            [t, indG, thetas, rhos, classes, iperf] = model_ThetaConded(pars(ipar,:), Gs, nT, nrr, trange, dt, rhobounds, mu_pri0);
            coefs(:,:,ipar) = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
            perfs(ipar) = iperf;

            fprintf('iIter=%3d/%3d | Theta Conded | ipar=%5d/%5d | iperf=%4.1f\n',iIter,nIter,ipar,npars,perfs(ipar))
        end

        dataStruct = struct("coefs",coefs,"pars",pars,"perfs",perfs);
        save(ifilename,'-fromstruct',dataStruct)

    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% RhoConded settings %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig_ms      = linspace(0.1,1,npar);
    sig_pri0s   = linspace(0.4,4,npar);
    sig_rs      = linspace(0.02,0.2,npar);
    dif_rs      = linspace(0.02,0.2,npar);
    
    pars        = allcomb(sig_ms,sig_pri0s,sig_rs,dif_rs);
    npars       = size(pars,1);

    isavename   = sprintf([dir0 '/result/models/RhoConded/'...
        'ms%.2f%.2f_pri0%.2f%.2f_srs%.2f%.2f_drs%.2f%.2f_nT%d_npar%d'],...
        sig_ms(1),sig_ms(end),sig_pri0s(1),sig_pri0s(end),sig_rs(1),sig_rs(end),dif_rs(1),dif_rs(end),nT,npar);
    if isempty(dir(isavename))
        mkdir(isavename)
    end
    ifilename = [isavename '/Iter' num2str(iIter) '.mat'];

    if isempty(dir(ifilename))

        coefs = NaN(8,3,npars); % ncoef, nbin, nvar, npar
        perfs = NaN(npars,1);
        for ipar = 1:npars
            [t, indG, thetas, rhos, classes, iperfs] = model_RhoConded(pars(ipar,:), Gs, nT, nrr, rrange0, rrange, rhobounds, mu_pri0);
            coefs(:,:,ipar) = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes);
            perfs(ipar) = iperfs;

            fprintf('iIter=%3d/%3d | RhoConded | ipar=%5d/%5d | iperf=%4.1f\n',iIter,nIter,ipar,npars,perfs(ipar))
        end

        dataStruct = struct("coefs",coefs,"pars",pars,"perfs",perfs);
        save(ifilename,'-fromstruct',dataStruct)

    end

end
