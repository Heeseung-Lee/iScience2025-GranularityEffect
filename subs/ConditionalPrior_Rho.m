function y = ConditionalPrior_Rho(iG,iclass,sigma,rangerho)

rhobound    = linspace(0,1,iG+1);
irange      = rhobound(iclass:iclass+1);
drho        = rangerho(2) - rangerho(1);
irho_uni    = rangerho>irange(1) & rangerho<irange(2);
iInd        = find(irho_uni==1);
nrho        = length(rangerho);

if sigma ~= 0
    irho        = zeros(1,nrho);
    nInd        = length(iInd);
    for i = 1:nInd
        kInd    = iInd(i);
        jrange  = rangerho - rangerho(kInd);
        irho    = irho + normpdf(jrange,0,sigma);
    end
else
    irho        = irho_uni;
end
irho    = irho/sum(irho*drho);
y       = range_truncate_from0to1(irho,rangerho);
y       = y/sum(y*drho);

end