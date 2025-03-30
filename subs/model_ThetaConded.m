function [t, indG, thetas, rhos, classes, perfs] = model_ThetaConded(pars, Gs, nT, nrr, trange, dt, rhobounds, mu_pri0)
    % Unpack parameters
    sig_m       = pars(1);
    sig_mm      = pars(2);
    sig_pri0    = pars(3);
    dif_t       = pars(4);

    % Initialize prior and convolution kernel
    ipri_t      = normpdf(trange, mu_pri0, sig_pri0); % initialize
    iconv       = normpdf(trange, 0, dif_t);

    % Generate random data
    t   = normrnd(0, 1, [nT, 1]);
    m   = normrnd(t, sig_m);
    mm  = normrnd(t, sqrt(sig_m^2 + sig_mm^2));

    % Initialize variables
    indG    = randi(2, [nT, 1]);
    thetas  = NaN(nT, 1);
    rhos    = NaN(nT, 1);
    classes = NaN(nT, 1);

    % Iterative computations
    for iT = 1:nT
        im          = m(iT);
        ilik_t      = normpdf(trange, im, sig_m);
        ilik_t(ilik_t == 0) = eps;
        ipri_t(ipri_t == 0) = eps;

        % Posterior computation
        imu_pos     = ilik_t .* ipri_t;
        imu_pos     = trange(find(imu_pos == max(imu_pos), 1)); % MAP
        cpri_t      = cumsum(ipri_t * dt);
        cpri_t(cpri_t > 1) = 1;
        ir_hat      = cpri_t(find(abs(trange - imu_pos) == min(abs(trange - imu_pos)), 1));

        % Save results
        thetas(iT)  = imu_pos;
        rhos(iT)    = ir_hat;

        % Determine class
        jG      = indG(iT);
        inG     = Gs(jG);
        jclass  = NaN;
        for ig = 1:inG
            irange = rhobounds{jG}(ig:ig+1);
            if ir_hat >= irange(1) && ir_hat <= irange(2)
                jclass      = ig;
                classes(iT) = jclass;
                break
            end
        end

        if isnan(jclass)
            pause
        end

        % Prior update
        im              = mm(iT);
        imu_pri         = (im * sig_pri0^2 + mu_pri0 * sig_mm^2) / (sig_pri0^2 + sig_mm^2);
        isig_pri        = sig_pri0 * sig_mm / sqrt(sig_pri0^2 + sig_mm^2);

        if iT ~= nT
            isd         = max([sig_m sig_mm isig_pri]);
            imax        = max(abs([im m(iT+1)]));
            trange      = (-2.5 * isd - imax + dt / 2):dt:(2.5 * isd + imax - dt / 2);
        end

        ipri_t          = normpdf(trange, imu_pri, isig_pri);
        cpri_t          = normcdf(trange, imu_pri, isig_pri);

        irrange         = rhobounds{jG}(jclass:jclass+1);
        iInd            = cpri_t > irrange(1) & cpri_t < irrange(2);
        ipri_t(~iInd)   = 0;
        ipri_t          = ipri_t / sum(ipri_t * dt);
        ipri_t          = conv(ipri_t, iconv, 'same');
        ipri_t          = ipri_t / sum(ipri_t * dt);

    end

    y0 = NaN(size(classes));
    y0(indG==1) = classes(indG==1)>1;
    y0(indG==2) = classes(indG==2)>4;
    perfs = mean(y0 == (t>0));
end
