function [jd1,jd2,jd3] = binarizer(d0,ig)

jd1 = zeros(size(d0));
jd2 = zeros(size(d0));
jd3 = zeros(size(d0));
gs = unique(ig(~isnan(ig(:))));
for i = gs
    jInd = ig==i;
    id0 = d0(jInd);

    id1 = zeros(size(id0));
    id2 = zeros(size(id0));
    id3 = zeros(size(id0));
    if i == 2
        id1(id0>=2)          = 1;
        id1(id0<2)           = 0;
    elseif i == 4
        id1(id0>=3)          = 1;
        id1(id0<3)           = 0;
        id2(id0==2|id0==4)    = 1;
        id2(id0==1|id0==3)    = 0;
    elseif i == 8
        id1(id0>=5)          = 1;
        id1(id0<5)           = 0;
        id2(id0==3|id0==4|id0==7|id0==8) = 1;
        id2(id0==1|id0==2|id0==5|id0==6) = 0;
        id3(id0==2|id0==4|id0==6|id0==8) = 1;
        id3(id0==1|id0==3|id0==5|id0==7) = 0;
    end
    id1(isnan(id0)) = NaN;
    id2(isnan(id0)) = NaN;
    id3(isnan(id0)) = NaN;

    jd1(jInd) = id1;
    jd2(jInd) = id2;
    jd3(jInd) = id3;
end