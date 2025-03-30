function [theta, theta_hat, rho_hat, classes] = model1_theta_only_update_3tasks(sig_m, sig_diff_theta, nIter, sig_theta, rangetheta, drho, rangerho, nrho, nT)

% stimulus sample is shared by the simulations
theta           = normrnd(0,sig_theta,[nT nIter]);
im              = normrnd(theta,sig_m);

% classification
theta_hat   = cell(1,3);
rho_hat     = cell(1,3);
classes     = cell(1,3);
for itask = 1:3
    switch itask
        case 1 % task1: G1 = 8, G0 = 2
            nGs     = [8 2];
        case 2 % task2: G1 = 2, G0 = 2
            nGs     = [2 2];
        case 3 % task3: G1 = 2, G0 = 8
            nGs     = [2 8];
    end
    rhobounds       = cell(nT,1);
    for iT = 1:nT
        iG              = nGs(iT);
        rhobounds{iT}   = linspace(0,1,iG+1);
    end
    
    imu0            = zeros(1,nIter);
    isig0           = 1;
    irho_pri        = ones(nIter,nrho)/drho/nrho;

    iclasses    = NaN(nT,nIter);
    ithetas     = NaN(nT,nIter);
    irho        = NaN(nT,nIter);
    for iT = 1:nT
        [itheta_hat, isig_pos]  = theta_inference(im(iT,:), imu0, isig0, sig_m);
        [iclass, irho_hat]      = rho_inference(nGs(iT), im(iT,:), sig_m, imu0, isig0, rangetheta, rangerho, drho, irho_pri, rhobounds{iT});
        [imu0, isig0]           = prior_update(itheta_hat, isig_pos, sig_diff_theta);

        iclasses(iT,:)          = iclass;
        ithetas(iT,:)           = itheta_hat;
        irho(iT,:)              = irho_hat;
    end
    theta_hat{itask}    = ithetas;
    rho_hat{itask}      = irho;
    classes{itask}      = iclasses;
end
end





