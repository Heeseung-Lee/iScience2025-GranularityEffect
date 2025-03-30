function ps = Convolve_Rho(ix,ix0,iy0,isigma)

nIter   = size(iy0,1);
ps      = NaN(size(iy0));
if isigma == 0
    ps  = iy0;
else
    for iIter = 1:nIter
        ix0         = round(ix0,5);
        ix          = round(ix,5);
        iInd        = find(ix0==ix(1)):find(ix0==ix(end));
        iy          = zeros(length(ix0),1);
        iy(iInd)    = iy0(iIter,:);
        dx          = ix0(2) - ix0(1);
        icumy       = cumsum(iy*dx);
        iInd        = find(icumy>0.0001,1,'first'):find(icumy<0.9999,1,'last');
        if isempty(iInd)
            iInd    = find(icumy>0.9999,1,'first');
        end
        nx          = length(ix0);
        ip          = zeros(1,nx);
        nInd        = length(iInd);
        for i = 1:nInd
            kInd    = iInd(i);
            jx      = ix0 - ix0(kInd);
            ip      = ip + iy(kInd)*normpdf(jx,0,isigma);
        end
        ip          = ip/sum(ip*dx);
        ip          = range_truncate_from0to1(ip,ix0);
        ps(iIter,:) = ip;
    end
end

end