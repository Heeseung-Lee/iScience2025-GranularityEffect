function [imu0, isig0] = prior_update(itheta_hat, isig_pos, sig_diff_theta)
imu0    = itheta_hat;
isig0   = sqrt(isig_pos^2 + sig_diff_theta^2);
end