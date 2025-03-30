function [y,thetabound] = ConditionalPrior_Theta(iG,iclass,sigma,rangetheta,range_theta)

thetabound  = linspace(range_theta(1),range_theta(2),iG+1);
irange      = thetabound(iclass:iclass+1);
dtheta      = rangetheta(2) - rangetheta(1);
itheta_uni  = rangetheta>=irange(1) & rangetheta<irange(2);
iInd        = find(itheta_uni==1);
ntheta      = length(rangetheta);

itheta      = zeros(1,ntheta);
nInd        = length(iInd);
for i = 1:nInd
    kInd    = iInd(i);
    jrange  = rangetheta - rangetheta(kInd);
    itheta  = itheta + normpdf(jrange,0,sigma);
end
y       = itheta/sum(itheta*dtheta);


% ntheta      = length(rangetheta);
% dtheta      = rangetheta(2) - rangetheta(1);
% iInd        = find(rangetheta>range_theta(1),1,'first'):find(rangetheta<range_theta(2),1,'last');
% itheta      = zeros(1,ntheta);
% nInd        = length(iInd);
% for i = 1:nInd
%     kInd    = iInd(i);
%     jrange  = rangetheta - rangetheta(kInd);
%     itheta  = itheta + normpdf(jrange,0,sigma);
% end
% jtheta      = itheta/sum(itheta*dtheta);
% icdf        = cumsum(jtheta*dtheta);
% nbound      = iG+1;
% pbound      = linspace(0.0001,0.9999,nbound);
% pInd        = NaN(1,nbound);
% for ibound  = 1:nbound
%     ie              = abs(icdf - pbound(ibound));
%     pInd(ibound)    = find(ie == min(ie),1);
% end
% 
% iInd        = pInd(iclass):pInd(iclass+1);
% itheta      = zeros(1,ntheta);
% nInd        = length(iInd);
% for i = 1:nInd
%     kInd    = iInd(i);
%     jrange  = rangetheta - rangetheta(kInd);
%     itheta  = itheta + normpdf(jrange,0,sigma);
% end
% y           = itheta/sum(itheta*dtheta);

end