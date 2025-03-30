function [t, indG, thetas, rhos, classes, perfs] = model_ThetaUpdate(pars, Gs, nT, rhobounds, mu_pri0)
    % Unpack parameters
    sig_m       = pars(1);
    sig_mm      = pars(2);
    sig_pri0    = pars(3);

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
        if iT == 1
            imu_pri  = mu_pri0;
            isig_pri = sig_pri0;
        end

        % Compute posterior mean and MAP estimate
        im      = m(iT);
        imu_pos = (im * isig_pri^2 + imu_pri * sig_m^2) / (isig_pri^2 + sig_m^2);
        ir_hat  = normcdf(imu_pos, imu_pri, isig_pri); % MAP

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

        % Prior update
        im        = mm(iT);
        imu_pri   = (im * sig_pri0^2 + mu_pri0 * sig_mm^2) / (sig_pri0^2 + sig_mm^2);
        isig_pri  = sig_pri0 * sig_mm / sqrt(sig_pri0^2 + sig_mm^2);
    end

    y0 = NaN(size(classes));
    y0(indG==1) = classes(indG==1)>1;
    y0(indG==2) = classes(indG==2)>4;
    perfs = mean(y0 == (t>0));
end
