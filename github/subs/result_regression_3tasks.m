function coef = result_regression_3tasks(theta, theta_hat, rho_hat, classes)

% task1: G1 = 8, G0 = 2
% task2: G1 = 2, G0 = 2
% task3: G1 = 2, G0 = 8

coef = NaN(7,3);
% 1: reproduction:      G1=8, G0=2
% 2: reproduction:      G1=2, G0=2
% 3: scaling:           G1=8, G0=2
% 4: scaling:           G1=2, G0=2
% 5: classification:    G1=8, G0=2
% 6: classification:    G1=2, G0=2
% 7: classification:    G1=2, G0=8

s1  = theta(1,:);
s0  = theta(2,:);

% reproduction
for i = 1:2
    x           = [classes{i}(1,:)' s0' s1'];
    y           = theta_hat{i}(2,:)';
    inan        = isnan(sum([y x],2));
    x           = zscore(x(~inan,:));
    y           = zscore(y(~inan));
    icoef       = glmfit(x,y);
    coef(i,:)   = icoef(2:end);
end

% scaling
for i = 1:2
    x           = [classes{i}(1,:)' s0' s1'];
    y           = rho_hat{i}(2,:)';
    y           = icdf('normal',y,0,2);
    inan        = isnan(sum([y x],2));
    x           = zscore(x(~inan,:));
    y           = zscore(y(~inan));
    icoef       = glmfit(x,y);
    coef(2+i,:) = icoef(2:end);
end

% classification
for i = 1:3
    x           = [classes{i}(1,:)' s0' s1'];
    y           = classes{i}(2,:)';
    inan        = isnan(sum([y x],2));
    x           = zscore(x(~inan,:));
    y           = zscore(y(~inan));
    icoef       = glmfit(x,y);
    coef(4+i,:) = icoef(2:end);
end

end