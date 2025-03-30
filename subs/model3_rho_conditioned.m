function results = model3_rho_conditioned(Gs, sig_m, sig_diff_theta, sig_diff_rho, nT, sig_theta, rangetheta, drho, rangerho, rangerho0, nrho)

nG          = length(Gs);
results     = cell(1,2);

for g = 1:nG

    iG  = Gs(g);
    
    % conditional prior of rho
    rhobound    = linspace(0,1,iG+1);
    prior_rho   = NaN(iG,nrho);
    for iclass = 1:iG
        prior_rho(iclass,:)     = ConditionalPrior_Rho(iG,iclass,sig_diff_rho,rangerho0);
    end

    irho_pri        = ones(1,nrho)/drho/nrho;
    imu0            = 0;
    isig0           = 1;
    iresults        = NaN(nT,5);
    for iT = 1:nT
        iclass = 0;
        while iclass == 0
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
            for ig = 1:iG
                irange          = rhobound(ig:ig+1);
                if irho_hat >= irange(1) && irho_hat < irange(2)
                    iclass  = ig;
                    break
                end
            end
        end

        % prior update
        imu0        = itheta_pos;
        isig0       = sqrt(isig_pos^2 + sig_diff_theta^2);
        irho_pri    = prior_rho(iclass,:);

        iresults(iT,:) = [itheta itheta_pos irho_likh irho_hat iclass];

    end

    results{g}  = iresults;
end

end