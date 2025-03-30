function ip = Convolve_Theta(ix,iy0,isigma)

if isigma == 0
    ip  = iy0;
else
    dx          = ix(2) -ix(1);
    icumy       = cumsum(iy0*dx);
    iInd        = find(icumy>0.0001,1,'first'):find(icumy<0.9999,1,'last');
    nx          = length(ix);
    ip          = zeros(1,nx);
    nInd        = length(iInd);
    for i = 1:nInd
        kInd    = iInd(i);
        jx      = ix - ix(kInd);
        ip      = ip + iy0(kInd)*normpdf(jx,0,isigma);
    end
    ip          = ip/sum(ip*dx);
end

end