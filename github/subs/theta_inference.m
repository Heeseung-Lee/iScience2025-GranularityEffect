function [itheta_hat, isig_pos] = theta_inference(im, imu0, isig0, sig_m)
itheta_hat  = (im*isig0^2 + imu0*sig_m^2)/(isig0^2 + sig_m^2);
isig_pos    = isig0*sig_m/sqrt(isig0^2 + sig_m^2);
end