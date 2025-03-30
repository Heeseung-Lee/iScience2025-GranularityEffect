function [coefs, warnings] = result_regression(results,nTB,imodel)

s0      = results(:,1); % physical size
t0      = results(:,2); % theta estimate
r0      = results(:,3); % rho estimate
c0      = results(:,4); % class estimate


s   = s0;
t   = t0;
r   = r0;
c   = c0;
name_s = cell(nTB,1);
name_t = cell(nTB,1);
name_c = cell(nTB,1);
for iTB = 1:nTB
    is  = [NaN(iTB,1); s0(1:end-iTB)];
    it  = [NaN(iTB,1); t0(1:end-iTB)];
    ir  = [NaN(iTB,1); r0(1:end-iTB)];
    ic  = [NaN(iTB,1); c0(1:end-iTB)];

    s   = cat(2,s,is);
    t   = cat(2,t,it);
    r   = cat(2,r,ir);
    c   = cat(2,c,ic);

    name_s{iTB} = ['stim_' num2str(iTB)];
    name_t{iTB} = ['repro_' num2str(iTB)];
    name_c{iTB} = ['class_' num2str(iTB)];
end
names  = [name_s; name_c];

ys          = [t0(:,1) r0(:,1) c0(:,1)];
ny          = size(ys,2);
coefs       = NaN(nTB*2,ny);
warnings    = ones(1,3);
if length(unique(results(:,4))) ~= 1
    for iy = 1:3
        y               = ys(:,iy);
        switch imodel
            case 1
                x       = [s(:,2:end) c(:,2:end)]; % x: stimulus and class
            case 2
                x       = [r(:,2:end) c(:,2:end)]; % x: reproduction and class
            case 3
                x       = [s(:,2:end) r(:,2:end) c(:,2:end)]; % x: stimulus, reproduction and class
        end
        reg             = [y x];
        reg             = zscore(reg(~isnan(sum(reg,2)),:));
        lastwarn('')
        coef            = glmfit(reg(:,2:end),reg(:,1));
        warnMsg         = lastwarn;
        if ~isempty(warnMsg)
            coef            = NaN(2*nTB+1,1);
            warnings(iy)    = 0;
        end
        coefs(:,iy)     = coef(2:end);
    end
end
coefs   = cell2table([names num2cell(coefs)],'VariableNames',{'reg','repro','scale','class'});
end