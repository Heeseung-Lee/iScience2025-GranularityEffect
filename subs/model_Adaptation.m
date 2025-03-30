function [t, indG, thetas, rhos, classes, perfs] = model_Adaptation(pars, Gs, nT, trange, dt, rhobounds, mu_pri0)
    % Unpack parameters
    sig_m       = pars(1);
    igain       = pars(2);
    sig_pri0    = pars(3);

    % Initialize prior and gain
    ipri_t      = normpdf(trange, mu_pri0, sig_pri0);
    ilik_gain   = ones(size(trange));

    % Generate random data
    t   = normrnd(0, 1, [nT, 1]);
    m   = normrnd(t, sig_m);

    % Initialize output variables
    indG    = randi(2, [nT, 1]);
    thetas  = NaN(nT, 1);
    rhos    = NaN(nT, 1);
    classes = NaN(nT, 1);

    % Iterative computations
    for iT = 1:nT
        im          = m(iT);
        ilik_t0     = normpdf(trange, im, sig_m);
        ilik_t      = ilik_t0 .* ilik_gain;

        % Posterior computation
        imu_pos     = ilik_t .* ipri_t;
        imu_pos     = trange(find(imu_pos == max(imu_pos), 1)); % MAP estimate

        % Cumulative prior update
        cpri_t      = cumsum(ipri_t * dt);
        cpri_t(cpri_t > 1) = 1;
        ir_hat      = cpri_t(find(abs(trange - imu_pos) == min(abs(trange - imu_pos)), 1)); % MAP

        % Save results
        thetas(iT)  = imu_pos;
        rhos(iT)    = ir_hat;

        % Determine class
        jG      = indG(iT);
        inG     = Gs(jG);
        jclass  = NaN;
        for ig = 1:inG
            irange = rhobounds{jG}(ig:ig+1);
            if ir_hat >= irange(1) && ir_hat < irange(2)
                jclass      = ig;
                classes(iT) = jclass;
                break;
            end
        end

        % Gain update
        ilik_gain   = 1 - igain * ilik_t0 / max(ilik_t0);
    end

    y0 = NaN(size(classes));
    y0(indG==1) = classes(indG==1)>1;
    y0(indG==2) = classes(indG==2)>4;
    perfs = mean(y0 == (t>0));
end
