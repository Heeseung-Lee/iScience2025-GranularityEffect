function [y,y1,x] = range_truncate_from0to1(y0,irange)

drho                = irange(2) - irange(1);
y0                  = y0/sum(y0*drho);
jInd                = [find(irange>0,1,'first') find(irange<1,1,'last')];
nrho                = length(irange);


while sum(y0(1:jInd(1)-1)) + sum(y0(jInd(2)+1:end)) ~= 0
    iy                  = zeros(1,nrho);
    iy(jInd(1) + (0:jInd(1)-2)) = fliplr(y0(1:jInd(1)-1));
    iy(jInd(1):nrho)    = y0(jInd(1):nrho) + iy(jInd(1):nrho);

    jy                  = zeros(1,nrho);
    jy(jInd(2) + (jInd(2)+1:nrho) - nrho) = fliplr(iy(jInd(2)+1:nrho));
    jy(1:jInd(2))       = iy(1:jInd(2)) + jy(1:jInd(2));

    y0 = jy;
end

y       = y0(jInd(1):jInd(2));
y       = y/sum(y*drho);
x       = irange(jInd(1):jInd(2));
y1      = y0(jInd(1):jInd(2));

end