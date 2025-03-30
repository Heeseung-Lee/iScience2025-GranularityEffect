function results = model2_theta_rho_update(Gs, sig_m, sig_diff_theta, ...
    sig_diff_rho, nT, sig_theta, rangetheta, drho, rangerho, rangerho0, nrho)

nG              = length(Gs);
rhobounds       = cell(nG,1);
for g = 1:nG
    iG              = Gs(g);
    rhobounds{g}    = linspace(0,1,iG+1);
end
results         = NaN(nT,4+nG);

irho_pri        = ones(1,nrho)/drho/nrho;
imu0            = 0;
isig0           = 1;
for iT = 1:nT
    iclasses    = zeros(1,nG);
    while sum(iclasses) == 0
        % theta generate
        itheta      = normrnd(0,sig_theta);
        im          = normrnd(itheta,sig_m);

        % theta inference
        itheta_pos  = (im*isig0^2 + imu0*sig_m^2)/(isig0^2 + sig_m^2);
        isig_pos    = isig0*sig_m/sqrt(isig0^2 + sig_m^2);

        % rho inference
        rhox        = normcdf(rangetheta,imu0,isig0);
        jrho_lik    = normpdf(rangetheta,im,sig_m);
        [~,iInd]    = unique(rhox);
        irho_lik    = interp1(rhox(iInd),jrho_lik(iInd),rangerho);
        irho_lik(isnan(irho_lik)) = 0;
        irho_lik    = irho_lik/sum(irho_lik*drho);
        irho_likh   = sum(rangerho.*irho_lik*drho);
        irho_pos    = irho_lik.*irho_pri;
        irho_pos    = irho_pos/sum(irho_pos*drho);
        irho_hat    = sum(rangerho.*irho_pos*drho);
        for g = 1:nG
            iG  = Gs(g);
            for ig = 1:iG
                irange          = rhobounds{g}(ig:ig+1);
                if irho_hat >= irange(1) && irho_hat < irange(2)
                    iclass      = ig;
                    iclasses(g) = iclass;
                    break
                end
            end
        end
    end

    % prior update
    imu0        = itheta_pos;
    isig0       = sqrt(isig_pos^2 + sig_diff_theta^2);
    irho_pri    = Convolve_Rho(rangerho,rangerho0,irho_pos,sig_diff_rho);

    results(iT,:)   = [itheta itheta_pos irho_likh irho_hat iclasses];

end

end