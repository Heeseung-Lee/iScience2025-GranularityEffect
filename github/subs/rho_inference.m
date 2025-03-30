function [iclass, irho_hat, irho_pos] = rho_inference(nG, im, sig_m, imu0, isig0, rangetheta, rangerho, drho, irho_pri, rhobounds)
nm          = length(im);
iclass      = NaN(nm,1);
irho_hat    = NaN(nm,1);
irho_pos    = NaN(nm,length(rangerho));

for i = 1:nm

    rhox            = normcdf(rangetheta,imu0(i),isig0);
    jrho_lik        = normpdf(rangetheta,im(i),sig_m);
    [~,iInd]        = unique(rhox);
    irho_lik        = interp1(rhox(iInd),jrho_lik(iInd),rangerho);
    irho_lik(isnan(irho_lik)) = 0;
    irho_lik        = irho_lik/sum(irho_lik*drho);
    irho_pos(i,:)   = irho_lik.*irho_pri(i,:);
    irho_pos(i,:)   = irho_pos(i,:)/sum(irho_pos(i,:)*drho);
    irho_hat(i)     = sum(rangerho.*irho_pos(i,:)*drho);

    for ig = 1:nG
        irange          = rhobounds(ig:ig+1);
        if irho_hat(i) >= irange(1) && irho_hat(i) < irange(2)
            iclass(i)   = ig;
            break
        end
    end
end
end