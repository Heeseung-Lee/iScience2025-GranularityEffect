

jresults        = load('/Volumes/CSNL_new/people/HSL/projects/granularity/results/Figures/humans/Boostrap/nIter10000/results_Scaling_ALL.mat');
iresults        = load('/Volumes/CSNL_new/people/HSL/projects/granularity/results/Figures/humans/Boostrap/nIter10000/results_PreCur.mat');

nIterFigure = 30;
nT          = 1000;
npar        = 10;
sig_ms      = linspace(0.1,1,npar);
sig_mms     = linspace(0.4,4,npar);
sig_pri0s   = linspace(0.4,4,npar);

isavename   = sprintf('/Volumes/CSNL_new/people/HSL/projects/granularity/results/Figures/models/BMBU_ThetaUpdate/ms%.2f%.2f_mms%.2f%.2f_pri0%.2f%.2f_nT%d_npar%d',...
    sig_ms(1),sig_ms(end),sig_mms(1),sig_mms(end),sig_pri0s(1),sig_pri0s(end),nT,npar);
acoefs = [];
for iIter = 1:nIterFigure
    if iIter == 1
        load([isavename '/Iter' num2str(iIter) '.mat'],'coefs','pars');
    else
        load([isavename '/Iter' num2str(iIter) '.mat'],'coefs');
    end
    acoefs = cat(4,acoefs,coefs);
end
coefs       = mean(acoefs,4,'omitnan');
ierror      = [abs(squeeze(coefs(:,1,:)) - jresults.coefs(1:6,2)); abs(squeeze(coefs(:,2,:)) - jresults.coefs(1:6,1)); abs(squeeze(coefs(:,3,:)) - iresults.coefs(1:6,1))];
ierror      = sum(ierror(1:18,:));
iInd_best   = find(ierror == min(ierror));
coef_best   = coefs(:,:,iInd_best);


