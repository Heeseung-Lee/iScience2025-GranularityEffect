function AP = get_AttributionProb(d0,s0,g0,nStimBin)

ng          = [2 8];
Binval      = linspace(0,1,nStimBin+1);
[~, rank]   = sort(s0, 1, 'ascend'); % 오름차순 정렬
[~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
s0          = order/size(s0,1);
pss         = NaN(size(s0));
for ibin = 1:nStimBin
    if ibin == nStimBin
        iInd    = s0>=Binval(ibin) & s0<=Binval(ibin+1);
    else
        iInd    = s0>=Binval(ibin) & s0<Binval(ibin+1);
    end
    pss(iInd)   = ibin;
end
s0          = pss;

d1 = [NaN; d0(1:end-1)];
s1 = [NaN; s0(1:end-1)];
g1 = [NaN; g0(1:end-1)];

AP = cell(2,2);
for ig = 1:2
    
    ing     = ng(ig);
    inTrans = NaN(ing,nStimBin,ing,nStimBin);
    for id1 = 1:ing
        for is1 = 1:nStimBin
            for is0 = 1:nStimBin
                for id0 = 1:ing
                    iInd = (d0==id0) & (d1==id1) & (s1==is1) & (s0==is0) & (g0==ig) & (g1==ig);
                    inTrans(id0,is0,id1,is1) = sum(iInd(:));
                end
            end
        end
    end

    ifroms0 = [];
    ifromd1 = [];
    ifroms1 = [];
    for id0 = 1:ing
        ipTrans = squeeze(inTrans(id0,:,:,:)); % s0,d1,s1

        ifroms0 = cat(1,ifroms0,mean(ipTrans./sum(ipTrans,1),[2 3],'omitnan')');
        ifromd1 = cat(1,ifromd1,mean(ipTrans./sum(ipTrans,2),[1 3],'omitnan'));
        ifroms1 = cat(1,ifroms1,squeeze(mean(ipTrans./sum(ipTrans,3),[1 2],'omitnan'))');
    end
    ifroms0(ifroms0==0) = eps;
    ifromd1(ifromd1==0) = eps;
    ifroms1(ifroms1==0) = eps;
    KL_s0 = sum(ifroms0 .* log(ifroms0 ./ (ones(size(ifroms0))./size(ifroms0,2))),2);
    KL_d10 = (sum(ifromd1 .* log(ifromd1 ./ (ones(size(ifromd1))./size(ifromd1,2))),2)./KL_s0)';
    KL_s10 = (sum(ifroms1 .* log(ifroms1 ./ (ones(size(ifroms1))./size(ifroms1,2))),2)./KL_s0)';
    
    AP{ig,1} = KL_d10;
    AP{ig,2} = KL_s10;
end

end