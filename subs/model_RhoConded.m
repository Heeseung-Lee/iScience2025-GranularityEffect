function [t, indG, thetas, rhos, classes, perfs] = model_RhoConded(pars, Gs, nT, nrr, rrange0, rrange, rhobounds, mu_pri0)

    % Unpack parameters
    sig_m       = pars(1);
    sig_pri0    = pars(2);
    sig_r       = pars(3);
    dif_r       = pars(4);
    
    % Generate random data
    t   = normrnd(0, 1, [nT, 1]);
    m   = normrnd(t, sig_m);

    % Compute prior_rho
    prior_rho = cell(1, max(Gs));
    for iG = Gs
        iprior_rho = NaN(iG, nrr);
        for iclass = 1:iG
            iprior_rho(iclass, :) = ConditionalPrior_Rho(iG, iclass, dif_r, rrange0);
        end
        prior_rho{iG} = iprior_rho;
    end
    
    % Initialize variables
    indG    = randi(2, [nT, 1]);
    ir_pri  = ones([1, nrr]);
    thetas  = NaN(nT, 1);
    rhos    = NaN(nT, 1);
    classes = NaN(nT, 1);

    % Iterative computations
    for iT = 1:nT
        if iT == 1
            imu_pri  = mu_pri0;
            isig_pri = sig_pri0;
        end

        % Compute posterior mean and variance of theta
        im      = m(iT);
        imu_pos = (im * isig_pri^2 + imu_pri * sig_m^2) / (isig_pri^2 + sig_m^2);
        ir      = normcdf(imu_pos, imu_pri, isig_pri);

        % Compute likelihood and posterior of rho
        ir_lik  = normpdf(rrange0, ir, sig_r);
        ir_lik  = range_truncate_from0to1(ir_lik, rrange0);
        ir_pos  = ir_lik .* ir_pri;
        % ir_hat  = rrange(find(ir_pos == max(ir_pos), 1)); % MAP estimate
        ir_pos  = ir_pos/sum(ir_pos)/(rrange(2)-rrange(1));
        ir_hat  = sum(ir_pos.*rrange*(rrange(2)-rrange(1))); % average following Wei & Stocker

        % Save results
        thetas(iT) = imu_pos;
        rhos(iT)   = ir_hat;

        % Determine class
        jG     = indG(iT);
        inG    = Gs(jG);
        jclass = NaN;
        for ig = 1:inG
            irange = rhobounds{jG}(ig:ig+1);
            if ir_hat >= irange(1) && ir_hat < irange(2)
                jclass      = ig;
                classes(iT) = jclass;
                break;
            end
        end

        % Update prior of rho
        ir_pri    = prior_rho{inG}(classes(iT), :);
    end

    y0 = NaN(size(classes));
    y0(indG==1) = classes(indG==1)>1;
    y0(indG==2) = classes(indG==2)>4;
    perfs = mean(y0 == (t>0));
end
