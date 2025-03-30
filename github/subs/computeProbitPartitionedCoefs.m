function [coefs,yhat,ig,results] = computeProbitPartitionedCoefs(t, indG, thetas, rhos, classes)

indG1   = [NaN; indG(1:end-1)];
d0      = [thetas rhos classes];
d1      = [NaN; classes(1:end-1)];
s0      = t;
s1      = [NaN; s0(1:end-1)];
Xs      = [];
d1s     = [];
bis     = [2 8];
for ig = 1:2
    if ig == 1
        iInd = indG == 1 & indG1 == 1;
    else
        iInd = indG == 1 & indG1 == 2;
    end
    %
    id1 = d1(iInd);
    [id01,id02,id03] = binarizer(id1,ones(size(id1))*bis(ig));

    X   = [d0(iInd,:) id01 s0(iInd) s1(iInd) id02 id03];
    
    X(:,3) = X(:,3)>1;
    if ig == 1
        X(:,[1 4:6]) = zscore_HL(X(:,[1 4:6]));
    else
        X(:,[1 4:8]) = zscore_HL(X(:,[1 4:8]));
    end
    Xs  = cat(1,Xs,[X ig*ones(size(X,1),1)]);
    d1s = cat(1,d1s,d1(iInd));
end
Xs(:,[1 4 5 6]) = zscore_HL(Xs(:,[1 4 5 6]));
id1 = Xs(:,7:8);
iInd = sum(id1,2)==0;
id1(~iInd,:) = zscore_HL(id1(~iInd,:));
Xs(:,7:8) = id1;

coefs = NaN(8,3);
yhat = NaN(size(Xs,1),3);
for ivar = 1:3
    X                   = [Xs(:,ivar) Xs(:,4:6) Xs(:,4:6).*(Xs(:,9)-1) Xs(:,7:8)];
    lastwarn('')
    switch ivar
        case 1 % reproduction
            icoef       = glmfit(X(:,2:end),X(:,1));
        case 2 % scaling
            icoef       = glmfit(X(:,2:end),X(:,1),'normal','link','probit');
        case 3 % classification
            icoef       = glmfit(X(:,2:end),X(:,1),'binomial','link','probit');
    end
    warnMsg             = lastwarn;
    if ~isempty(warnMsg)
        icoef           = NaN(1,9);
    else
        switch ivar
            case 1
                yhat(:,ivar)    = glmval(icoef,X(:,2:end),'identity');
            otherwise
                yhat(:,ivar)    = glmval(icoef,X(:,2:end),'probit');
        end
    end
    coefs(:,ivar)   = icoef(2:end);
end
ig = Xs(:,end);

results.g       = ig;
results.y       = Xs(:,1:3);
results.yhat    = yhat;
results.d1      = d1s;
results.s0      = Xs(:,5);
results.s1      = Xs(:,6);

end