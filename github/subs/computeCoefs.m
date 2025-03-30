function coefs = computeCoefs(t, indG, thetas, rhos, classes, Gs, nStimBins, nBins)
    % computeCoefs: Calculate coefficients using logistic regression
    %
    % Inputs:
    % - t: time series data (nT x 1 vector)
    % - indG: group indices (nT x 1 vector)
    % - thetas, rhos, classes: feature variables (nT x 1 vectors)
    % - Gs: array of group sizes
    % - nStimBins: array of stimulus bins
    % - nBins: number of bins
    % - ipar: parameter index
    % - d0, d1: auxiliary matrices for logistic regression
    %
    % Output:
    % - coefs: calculated coefficients (6 x nBins x 2 x 1 array)

    indG1   = [NaN; indG(1:end-1)];
    d0      = [thetas>0 rhos>0.5 classes>Gs(indG)'/2];
    d1      = [NaN; classes(1:end-1)];

    % Initialize variables
    coefs = NaN(6, nBins, 2); % Adjust size as needed for your requirements
    
    for iBin = 1:nBins
        nStimBin = nStimBins(iBin);
        Binval = linspace(0, 1, nStimBin+1);

        % Normalize and bin the time series
        s0 = t;
        [~, rank] = sort(s0, 1, 'ascend');
        [~, order] = sort(rank, 1, 'ascend');
        s0 = order / size(s0, 1);
        pss = NaN(size(s0));
        for ibin = 1:nStimBin
            if ibin == nStimBin
                iInd = s0 >= Binval(ibin) & s0 <= Binval(ibin+1);
            else
                iInd = s0 >= Binval(ibin) & s0 < Binval(ibin+1);
            end
            pss(iInd) = ibin;
        end
        s0 = pss;
        stimvals = unique(s0(:));
        s1 = [NaN; s0(1:end-1)];

        % Initialize variables for groups
        ys = cell(1, 2);
        nn = cell(1, 2);
        xs = cell(1, 2);
        gs = [];

        for ig = 1:2
            if ig == 1
                iInd = indG == 1 & indG1 == 1;
            else
                iInd = indG == 1 & indG1 == 2;
            end

            id0 = d0(iInd, :);
            id1 = d1(iInd);
            is1 = s1(iInd);

            ni = nStimBin;
            nj = Gs(ig);

            kd0 = NaN(ni, nj, 3);
            nd0 = NaN(ni, nj);
            for i = 1:ni
                for j = 1:nj
                    iInd = find((is1 == stimvals(i)) & (id1 == j) & ~isnan(is1 + id1 + sum(id0, 2)));
                    jd0 = id0(iInd, :);
                    kd0(i, j, :) = sum(jd0, 1);
                    nd0(i, j) = length(iInd);
                end
            end
            ip = kd0 ./ nd0;
            in = nd0;

            ing = size(ip, 2);

            y = [mean(ip(:, 1:floor(ing/2), :), 2), mean(ip(:, ceil(1+ing/2):end, :), 2)];
            in = [sum(in(:, 1:floor(ing/2)), 2), sum(in(:, ceil(1+ing/2):end), 2)];
            x1 = repmat((1:nStimBin)', 1, 2);
            x2 = repmat(1:2, nStimBin, 1);
            nn{ig} = in(:);

            ys{ig} = reshape(y, size(y, 1) * size(y, 2), size(y, 3));
            xs{ig} = [x1(:), x2(:)];
            gs = cat(1, gs, ig * ones(size(y, 1) * size(y, 2), 1));
        end

        inan = isnan(sum(ys{1} + ys{2}, 2));
        nn{1}(inan) = NaN;
        nn{2}(inan) = NaN;
        in = [nn{1}; nn{2}];
        iy = [ys{1}; ys{2}];
        jx = [ones(size([xs{1}; xs{2}], 1), 1), [xs{1}; xs{2}]];
        ix = [jx, jx .* (gs - 1)];

        iX = [];
        for it = 1:length(in)
            jn = in(it);
            if jn > 0
                iX = cat(1, iX, repmat(ix(it, :), jn, 1));
            end
        end

        for ivar = 1:2
            iY = [];
            for it = 1:length(in)
                jn = in(it);
                if jn > 0
                    jp = iy(it, ivar);
                    jy = zeros(jn, 1);
                    jy(1:round(jn * jp)) = 1;
                    iY = cat(1, iY, jy);
                end
            end
            lastwarn('');
            if ~isempty(iX)
                r = glmfit(iX, iY, 'binomial', 'link', 'logit', 'constant', 'off', 'Options', statset('Display', 'final', 'MaxIter', 5));
            else
                r = NaN;
            end
            warnMsg = lastwarn;
            if ~isempty(warnMsg)
                r = NaN;
            end
            coefs(:, iBin, ivar) = r;
        end
    end
end
