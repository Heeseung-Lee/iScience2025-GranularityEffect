function rs = computeBetas(xs, stimvals, xvals, ng, nStimBins, iFigureCond, iBoot, ngs, ni, nj, iBin)
    % Input Arguments:
    % xs: Cell array containing input data for each group
    % stimvals: Vector of stimulus values
    % xvals: Vector of x values for z-scoring
    % ng: Vector specifying the number of conditions per group
    % nStimBins: Vector specifying the number of stimulus bins
    % iFigureCond: Boolean for whether to generate plots
    % iBoot: Integer specifying bootstrapping iteration
    % ngs: Total number of groups
    % ni, nj, nk: Dimensions for loop iterations
    % iBin, ifeature: Indices for organizing output

    betas = cell(ngs, 3);
    rs = NaN(3,ngs+1);
    xvals = zscore(xvals);
    
    for ig = 1:ngs
        nT = size(xs{ig}, 1);
        
        if iBoot == 1
            iInd0 = 1:nT;
        else
            iInd0 = randi(nT, [nT, 1]);
        end
        
        nk = ng(ig);
        kd0 = NaN(ni, nj, nk);
        nd0 = NaN(ni, nj, nk);
        
        for i = 1:ni
            for j = 1:nj
                for k = 1:nk
                    iy0 = xs{ig}(iInd0, 1);
                    is0 = xs{ig}(iInd0, 2);
                    is1 = xs{ig}(iInd0, 3);
                    id1 = xs{ig}(iInd0, 4);

                    iInd = find((is0 == stimvals(i)) & (is1 == stimvals(j)) & (id1 == k) & ~isnan(is0 + is1 + id1 + iy0));
                    kd0(i, j, k) = sum(iy0(iInd));
                    nd0(i, j, k) = length(iInd);
                end
            end
        end
        
        ip = kd0 ./ nd0;
        in = nd0;
        nStimBin = nStimBins(iBin);
        ing = size(ip, 3);
        
        % Compute s0 effect
        for is1 = 1:nStimBin
            for ic1 = 1:ing
                [betas{ig, 3}(is1, ic1)] = computeEffect(ip(:, is1, ic1), in(:, is1, ic1), xvals);
            end
        end

        % Compute s1 effect and c1 effect
        for is0 = 1:nStimBin
            if iFigureCond && iBoot == 1
                subplot(ngs, 1 + nStimBin, (1 + nStimBin) * (ig - 1) + is0);
                hold on;
                cl = jet(ing);
            end
            
            for i = 1:ing
                [betas{ig, 2}(i, is0), iconst] = computeEffect(ip(is0, :, i), in(is0, :, i), xvals);
                
                if iFigureCond && iBoot == 1
                    plotEffect(ip(is0, :, i), [iconst betas{ig, 2}(i, is0)],xvals, cl(i, :));
                end
            end
            
            if iFigureCond && iBoot == 1
                ylim([0, 1]);
                grid on;
            end
            
            for i = 1:nStimBin
                [betas{ig, 1}(i, is0)] = computeEffect(ip(is0, i, :), in(is0, i, :), zscore(1:ing));
            end
        end
    end

    % Summarize results
    for ic1s1s0 = 1:3
        iy = [];
        ix = [];
        ys = [];
        for ig = 1:ngs
            jy = betas{ig, ic1s1s0}(:);
            ys = cat(2, ys, mean(jy, 'omitnan'));
            iy = cat(1, iy, jy);
            ix = cat(1, ix, ng(ig) * ones(size(jy)));
        end
        ir = glmfit(ix, iy);
        rs(ic1s1s0, :) = [ys, ir(2)];
    end
end

function [beta, const] = computeEffect(ip, in, xvals)
    % Helper function to compute regression coefficients
    iX = [];
    iY = [];
    for it = 1:length(in)
        jn = in(it);
        if jn > 0
            jp = ip(it);
            jy = zeros(jn, 1);
            jy(1:round(jn * jp)) = 1;
            iY = cat(1, iY, jy);
            iX = cat(1, iX, repmat(xvals(it), jn, 1));
        end
    end
    try
        r = glmfit(iX, iY, 'binomial', 'link', 'logit', 'LikelihoodPenalty', 'jeffreys-prior');
        beta = r(2);
        const = r(1);
    catch
        beta = NaN;
        const = NaN;
    end
end

function plotEffect(ip, icoef, xvals, color)
    % Helper function to plot effects
    plot(xvals, ip, '-', 'color', color, 'linewidth', 1);
    yhat = glmval(icoef', xvals, 'logit');
    plot(xvals, yhat, '--o', 'color', color, 'linewidth', 1);
end
