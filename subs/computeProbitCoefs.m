function [coefs,yhat,ig,results] = computeProbitCoefs(t, indG, thetas, rhos, classes)

indG1   = [NaN; indG(1:end-1)];
d0      = [thetas rhos classes];
d1      = [NaN; classes(1:end-1)];
s0      = t;
s1      = [NaN; s0(1:end-1)];
Xs      = [];
d1s     = [];
for ig = 1:2
    if ig == 1
        iInd = indG == 1 & indG1 == 1;
    else
        iInd = indG == 1 & indG1 == 2;
    end
    %
    X   = [d0(iInd,:) d1(iInd) s0(iInd) s1(iInd)];
    X(:,3) = X(:,3)>1;
    X(:,[1 4 5 6]) = zscore_HL(X(:,[1 4 5 6]));
    Xs  = cat(1,Xs,[X ig*ones(size(X,1),1)]);
    d1s = cat(1,d1s,d1(iInd));
end
Xs(:,[1 4 5 6]) = zscore_HL(Xs(:,[1 4 5 6]));
coefs = NaN(6,3);
yhat = NaN(size(Xs,1),3);
for ivar = 1:3
    X                   = [Xs(:,ivar) Xs(:,4:6) Xs(:,4:6).*(Xs(:,7)-1)];
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
        icoef           = NaN(1,7);
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