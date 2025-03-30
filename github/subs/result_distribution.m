function [dist_t, dist_r, dist_c] = result_distribution(results, nbin, xbin, rangerho, nrho, inG, ntbin, tbin, nT)

s0          = results(:,1);
t0          = results(:,2);
r0          = results(:,3);
d0          = results(:,4);

i0          = NaN(nT,1);
for ibin = 1:nbin
    iInd        = s0>=xbin(ibin) & s0<xbin(1+ibin);
    i0(iInd)    = ibin;
end
i1          = [NaN; i0(1:end-1)];
d1          = [NaN; d0(1:end-1)];

dist_t  = NaN(inG,ntbin,nbin);
dist_r  = NaN(inG,nrho,nbin);
dist_c  = NaN(inG,inG,nbin);
for is1 = 1:nbin
    for id1 = 1:inG
        iInd                = i1==is1 & d1 == id1;

        it0                 = t0(iInd);
        ih                  = hist(it0,tbin);
        dist_t(id1,:,is1)   = ih;

        ir0                 = r0(iInd);
        ih                  = hist(ir0,rangerho);
        dist_r(id1,:,is1)   = ih;

        id0                 = d0(iInd);
        ih                  = hist(id0,1:inG);
        dist_c(id1,:,is1)   = ih;
    end
end
for is1 = 1:nbin
    idist               = squeeze(dist_t(:,:,is1));
    idist               = idist./sum(idist*(tbin(2)-tbin(1)),2,'omitnan');
    dist_t(:,:,is1)     = idist;

    idist               = squeeze(dist_r(:,:,is1));
    idist               = idist./sum(idist*(rangerho(2)-rangerho(1)),2,'omitnan');
    dist_r(:,:,is1)     = idist;

    idist               = squeeze(dist_c(:,:,is1));
    idist               = idist./sum(idist,2,'omitnan');
    dist_c(:,:,is1)     = idist;
end
dist_t  = mean(dist_t,3,'omitnan');
dist_r  = mean(dist_r,3,'omitnan');
dist_c  = mean(dist_c,3,'omitnan');

end