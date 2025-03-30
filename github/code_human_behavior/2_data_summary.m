clear all; clc; close all

downloadfolder = '/Users/heeseunglee/Documents/project/granularity/github';
savedir     = [downloadfolder '/result/humans'];
if isempty(dir(savedir))
    mkdir(savedir)
end
cd(savedir)

addpath([downloadfolder '/subs'])
addpath([downloadfolder '/data'])

nTBs            = [10 7];
tt              = {'Ring','Pitch'};
results         = cell(9,1);
results{1}      = load(['results_RingPitch_nTB' num2str(nTBs(1)) '.mat']);
results{2}      = load('results_GeneralMagnitude.mat');
results{3}      = load('results_PreCur.mat');
results{4}      = load('results_Scaling_LongStimNoMask.mat');
results{5}      = load('results_Scaling.mat');
results{6}      = load('results_Motor.mat');
results{7}      = load('results_ITI.mat');
results{8}      = load(['results_RingPitch_nTB' num2str(nTBs(2)) '_Individual.mat']);
results{10}     = load("data_ringpitch_seperate.mat");

fz          = 15;
lw          = 1.2;
dx          = [-1.5 -0.5 0.5 1.5]*0.05;
colors      = lines(6);
pCI         = 95;

pval_fdr    = cell(8,1);
r_ci        = cell(8,1);
r           = cell(8,1);
for iexp = 1:8
    if iexp == 7
        r{iexp}         = results{iexp}.acoefs;
        r_ci{iexp}      = results{iexp}.ar_ci;
        pval_fdr{iexp}  = results{iexp}.apval_fdr;
    else
        r{iexp}         = results{iexp}.coefs;
        r_ci{iexp}      = results{iexp}.r_ci;
        pval_fdr{iexp}  = results{iexp}.pval_fdr;
    end
end
%
n = 256; % colormap의 단계 수
top = [1 1 0]; % 노란색
middle = [0.5 0.5 0.5]; % 회색
bottom = [0 0 0]; % 파란색

% 상단부터 중간까지의 색상 전환
top_half = [linspace(bottom(1), middle(1), n/2)', ...
            linspace(bottom(2), middle(2), n/2)', ...
            linspace(bottom(3), middle(3), n/2)'];

% 중간부터 하단까지의 색상 전환
bottom_half = [linspace(middle(1), top(1), n/2)', ...
               linspace(middle(2), top(2), n/2)', ...
               linspace(middle(3), top(3), n/2)'];

% colormap 합치기
custom_colormap = [top_half; bottom_half];

% colormap 적용
colormap(custom_colormap);

%%


% Exp 1,2 - Ring, Pitch average

figure(1)
clf
nTB = nTBs(1);
for ifeature = 1:2
    subplot(6,4,ifeature)
    hold on
    for icond = [3 4 1 2]
        if icond == 4
            iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB)-1);
            x       = [1:nTB nTB+2:(1+2*nTB)];
        else
            iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB));
            x       = (1:(1+2*nTB));
        end
        icoefs  = squeeze(r{1}(iInd,ifeature));
        ici     = squeeze(r_ci{1}(iInd,ifeature,:));
        ipval   = squeeze(pval_fdr{1}(iInd,ifeature));
        for ixrange = 1:3
            switch ixrange
                case 1
                    xInd = find(x<=nTB);
                case 2
                    xInd = find(x==nTB+1);
                case 3
                    xInd = find(x>nTB+1);
            end
            incoef  = length(xInd);
            if max(icoefs(xInd)) > 1.4
                dy = -0.8;
            elseif max(icoefs(xInd)) > 0.6
                dy = -0.6;
            else
                dy = 0;
            end
            for ix = 1:incoef
                plot([x(xInd(ix)) x(xInd(ix))]+dx(icond),ici(xInd(ix),:)+dy,'_-','color',colors(icond,:),'linewidth',lw)
            end
            plot(x(xInd)+dx(icond),icoefs(xInd)+dy,'o-','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
            ip = ipval(xInd)<0.05;
            plot(x(xInd(ip))+dx(icond),icoefs(xInd(ip))+dy,'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
        end
    end
    plot([0 2*nTB+2],[0 0],'-k')
    xlim([0 2*nTB+2])
    ylim([-0.6 1])
    title(tt{ifeature})
    if nTB == 10
        set(gca,'xtick',[1 5 10 11 12 16 21],'xticklabel',[1 5 10 0 1 5 10],'ytick',[-0.3:0.3:0.3 0.6 0.9],'yticklabel',[-0.3:0.3:0.3 1.2 1.7],'fontsize',fz)
    else
        set(gca,'xtick',[1 nTB nTB+1 nTB+2 2*nTB+1],'xticklabel',[1 nTB 0 1 nTB],'ytick',[-0.3:0.3:0.3 0.6 0.9],'yticklabel',[-0.3:0.3:0.3 1.2 1.7],'fontsize',fz)
    end
end

% psychometric curve

nBin    = 5;
xval    = linspace(0,1,nBin+1);
ngs     = [2 4 8];
r2 = NaN(1,2);
for ifeature = 1:2
    for ivar = 1:2
        for ig = 1:3
            subplot(18,16,48+3*(ifeature-1)+16*(ivar-1)+ig)
            cla;
            hold on
            iInd    = results{1}.xy{ifeature}.g==ig;
            iy      = results{1}.xy{ifeature}.y(iInd);
            iyhat   = results{1}.xy{ifeature}.yhat(iInd);
            switch ivar
                case 1
                    ix = results{1}.xy{ifeature}.s0(iInd);
                case 2
                    ix = results{1}.xy{ifeature}.s1(iInd);
            end
            ic1     = results{1}.xy{ifeature}.d1(iInd);

            [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
            [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
            ip0         = order/size(ix,1);
            jx = NaN(size(ix));
            for iBin = 1:nBin
                if iBin == nBin
                    iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                else
                    iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                end
                jx(iInd) = iBin;
            end

            ing = ngs(ig);
            cl = jet(ing);
            for jc1 = 1:ing
                ky      = NaN(nBin,1);
                kyhat   = NaN(nBin,1);
                for iBin = 1:nBin
                    iInd        = (jx==iBin) & (ic1==jc1);
                    jy          = iy(iInd);
                    ky(iBin)    = mean(jy);
                    kyhat(iBin) = mean(iyhat(iInd));
                end
                plot(kyhat,'-','color',cl(jc1,:),'linewidth',1.2)
                plot(ky,'o','color','k','markerfacecolor',cl(jc1,:))

                set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
            end
            ylim([0 1])
            xlim([0.5 nBin+0.5])
        end
    end
    r2(ifeature) = corr(results{1}.xy{ifeature}.y,results{1}.xy{ifeature}.yhat)^2;
end

% Exp 1,2 - Ring, Pitch individual

nTB = nTBs(2);
for ifeature = 1:2
    subplot(6,4,2+ifeature)
    cla;
    hold on
    for icond = [3 4 1 2]
        if icond == 4
            iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB)-1);
            x       = [1:nTB nTB+2:(1+2*nTB)];
        else
            iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB));
            x       = (1:(1+2*nTB));
        end
        icoefs  = squeeze(r{8}(iInd,ifeature));
        ici     = squeeze(r_ci{8}(iInd,ifeature,:));
        ipval   = squeeze(pval_fdr{8}(iInd,ifeature));
        for ixrange = 1:3
            switch ixrange
                case 1
                    xInd = find(x<=nTB);
                case 2
                    xInd = find(x==nTB+1);
                case 3
                    xInd = find(x>nTB+1);
            end
            incoef  = length(xInd);
            if max(icoefs(xInd)) > 2
                dy = -1.5;
            elseif max(icoefs(xInd)) > 0.6
                dy = -1;
            else
                dy = 0;
            end
            for ix = 1:incoef
                plot([x(xInd(ix)) x(xInd(ix))]+dx(icond),ici(xInd(ix),:)+dy,'_-','color',colors(icond,:),'linewidth',lw)
            end
            plot(x(xInd)+dx(icond),icoefs(xInd)+dy,'o-','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
            ip = ipval(xInd)<0.05;
            plot(x(xInd(ip))+dx(icond),icoefs(xInd(ip))+dy,'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
        end
    end
    plot([0 2*nTB+2],[0 0],'-k')
    xlim([0 2*nTB+2])
    ylim([-0.9 1.2])
    title(tt{ifeature})
    if nTB == 10
        set(gca,'xtick',[1 5 10 11 12 16 21],'xticklabel',[1 5 10 0 1 5 10],'ytick',[-0.9:0.3:0.3 0.6 0.9],'yticklabel',[-0.9:0.3:0.3 1.6 2.4],'fontsize',fz)
    else
        set(gca,'xtick',[1 nTB nTB+1 nTB+2 2*nTB+1],'xticklabel',[1 nTB 0 1 nTB],'ytick',[-0.9:0.3:0.3 0.6 0.9],'yticklabel',[-0.9:0.3:0.3 1.6 2.4],'fontsize',fz)
    end
end

%%

% psychometric curve

nBin    = 3;
xval    = linspace(0,1,nBin+1);
ngs     = [2 4 8];
iSub    = 51; % representative subject: 19; 25; 51*
r2s = NaN(58,2);
for ifeature = 1:2
    for iSub = 1:58
        iy      = results{8}.xy{iSub,ifeature}.y;
        iyhat   = results{8}.xy{iSub,ifeature}.yhat;
        r2s(iSub,ifeature) = corr(iy,iyhat)^2;
    end
    for ivar = 1:2
        for ig = 1:3
            subplot(18,16,56+3*(ifeature-1)+16*(ivar-1)+ig)
            cla;
            hold on
            iInd    = results{8}.xy{iSub,ifeature}.g==ig;
            iy      = results{8}.xy{iSub,ifeature}.y(iInd);
            iyhat   = results{8}.xy{iSub,ifeature}.yhat(iInd);
            switch ivar
                case 1
                    ix = results{8}.xy{iSub,ifeature}.s0(iInd);
                case 2
                    ix = results{8}.xy{iSub,ifeature}.s1(iInd);
            end
            ic1     = results{8}.xy{iSub,ifeature}.d1(iInd);

            [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
            [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
            ip0         = order/size(ix,1);
            jx = NaN(size(ix));
            for iBin = 1:nBin
                if iBin == nBin
                    iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                else
                    iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                end
                jx(iInd) = iBin;
            end

            ing = ngs(ig);
            cl = jet(ing);
            for jc1 = 1:ing
                ky      = NaN(nBin,1);
                kyhat   = NaN(nBin,1);
                for iBin = 1:nBin
                    iInd        = (jx==iBin) & (ic1==jc1);
                    jy          = iy(iInd);
                    if ~isempty(jy)
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                end
                ix = find(~isnan(ky));
                plot(ix,kyhat(ix),'-','color',cl(jc1,:),'linewidth',1.2)
                plot(ix,ky(ix),'o','color','k','markerfacecolor',cl(jc1,:))

                set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
            end
            ylim([0 1])
            xlim([0.5 nBin+0.5])
        end
    end
end

[~,~,ci] = ttest(r2s);
mean(r2s)
ci(2,:)-mean(r2s)
%% bar graph

for ifeature = 1:2
    ir = [];
    for icond = [1 2]
        iInd    = (1+(icond-1)*(1+2*nTB)):(icond*(1+2*nTB));
        icoefs  = squeeze(results{8}.irs(iInd,ifeature,:))';
        ir = cat(2,ir,icoefs(:,[1 nTB+2]));
    end
    ir = [ir(:,1) ir(:,3) ir(:,2) ir(:,4)];
    im = mean(ir);
    [~,~,ici] = ttest(ir);

    subplot(8,10,20 + 5*ifeature)
    cla;
    hold on
    bar(im,'FaceColor','w','EdgeColor','k','linewidth',lw)
    for i = 1:4
        plot([i i]+normrnd(0,0.05,[size(ir,1),1]),ir(:,i),'.','color',[1 1 1]*0.6)
        plot([i i],ici(:,i),'k_-','linewidth',1.5,'markersize',11)
    end
    set(gca,'xtick',1:4,'xticklabel',{'c1','c1xg','s1','s1xg'},'ytick',-1:0.5:1,'fontsize',fz)
    if ifeature == 1
        ylim([-1 1])
    else
        ylim([-1.5 1.5])
    end
end

%%

% Exp 3 - Generel magnitude

dx = [-0.5 0.5]*0.2;
c = 1;
yt = {-0.8:0.4:0.8,[],[],[],[],[]};
tt = {'R2R','P2P','P2R','R2P','Same','Diff'};
for iseq = [1 4 2 3 5 6]
    subplot(9,15,45+c)
    hold on
    for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
        iInd    = [1+(icond-1)*3 icond*3];
        x       = [1 2];
        icoefs  = squeeze(r{2}(iInd,iseq));
        ici     = squeeze(r_ci{2}(iInd,iseq,:));
        ipval   = squeeze(pval_fdr{2}(iInd,iseq));
        for ixrange = 1:length(x)
            plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:),'linewidth',lw)
        end
        plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
        ip = ipval<0.05;
        plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
    end
    plot([0.5 3.5],[0 0],'-k')
    xlim([0.5 2.5])
    ylim([-1 1])
    title(tt{c})
    set(gca,'fontsize',fz,'xtick',1:2,'ytick',yt{c},'xticklabel',[])
    c = c + 1;
end

% psychometric curve

nBin    = 3;
xval    = linspace(0,1,nBin+1);
ngs     = [2 4];
iInds   = [1 3 4 2 5 6];
c = 1;
for ifeature = 1:3
    for pfeature = 1:2
        for ivar = 1:2
            for ig = 1:2
                subplot(18,16,128+2*(iInds(c)-1)+16*(ivar-1)+ig)
                cla;
                hold on
                iInd    = results{2}.xy{ifeature,pfeature}.g==ig;
                iy      = results{2}.xy{ifeature,pfeature}.y(iInd);
                iyhat   = results{2}.xy{ifeature,pfeature}.yhat(iInd);
                switch ivar
                    case 1
                        ix = results{2}.xy{ifeature,pfeature}.s0(iInd);
                    case 2
                        ix = results{2}.xy{ifeature,pfeature}.s1(iInd);
                end
                ic1     = results{2}.xy{ifeature,pfeature}.d1(iInd);

                [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
                [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
                ip0         = order/size(ix,1);
                jx = NaN(size(ix));
                for iBin = 1:nBin
                    if iBin == nBin
                        iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                    else
                        iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                    end
                    jx(iInd) = iBin;
                end

                ing = ngs(ig);
                cl = jet(ing);
                for jc1 = 1:ing
                    ky      = NaN(nBin,1);
                    kyhat   = NaN(nBin,1);
                    for iBin = 1:nBin
                        iInd        = (jx==iBin) & (ic1==jc1);
                        jy          = iy(iInd);
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                    plot(kyhat,'-','color',cl(jc1,:),'linewidth',1.2)
                    plot(ky,'o','color','k','markerfacecolor',cl(jc1,:))

                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                end
                ylim([0 1])
                xlim([0.5 nBin+0.5])
            end
        end
        c = c + 1;
    end
end

% Exp 4 - PreCur

tt = {'G1 vary','G0 vary'};
c = 1;
yt = {-0.3:0.3:0.3,[],[],[]};
for iseq = 1:2
    subplot(9,15,56+c)
    hold on
    for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
        iInd    = [1+(icond-1)*3 icond*3];
        x       = 1:2;
        icoefs  = squeeze(r{3}(iInd,iseq));
        ici     = squeeze(r_ci{3}(iInd,iseq,:));
        ipval   = squeeze(pval_fdr{3}(iInd,iseq));
        for ixrange = 1:length(x)
            plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:),'linewidth',lw)
        end
        plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
        ip = ipval<0.05;
        plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
    end
    plot([0.5 2.5],[0 0],'-k')
    xlim([0.5 2.5])
    ylim([-0.5 0.4])
    title(tt{iseq})
    set(gca,'fontsize',fz,'xtick',1:2,'ytick',yt{c},'xticklabel',[])
    c = c + 1;
end

% psychometric curve

nBin    = 4;
xval    = linspace(0,1,nBin+1);
for iseq = 1:2
    for ivar = 1:2
        for ig = 1:2
            subplot(18,16,138+2*iseq+16*(ivar-1)+ig)
            cla;
            hold on
            iInd    = results{3}.xy{iseq}.g==ig;
            iy      = results{3}.xy{iseq}.y(iInd);
            iyhat   = results{3}.xy{iseq}.yhat(iInd);
            switch ivar
                case 1
                    ix = results{3}.xy{iseq}.s0(iInd);
                case 2
                    ix = results{3}.xy{iseq}.s1(iInd);
            end
            ic1     = results{3}.xy{iseq}.d1(iInd);

            [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
            [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
            ip0         = order/size(ix,1);
            jx = NaN(size(ix));
            for iBin = 1:nBin
                if iBin == nBin
                    iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                else
                    iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                end
                jx(iInd) = iBin;
            end
            
            ing = length(unique(ic1));
            cl = jet(ing);
            for jc1 = 1:ing
                ky      = NaN(nBin,1);
                kyhat   = NaN(nBin,1);
                for iBin = 1:nBin
                    iInd        = (jx==iBin) & (ic1==jc1);
                    jy          = iy(iInd);
                    if ~isempty(jy)
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                end
                ix = find(~isnan(ky));
                plot(ix,kyhat(ix),'-','color',cl(jc1,:),'linewidth',1.2)
                plot(ix,ky(ix),'o','color','k','markerfacecolor',cl(jc1,:))
                set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
            end
            ylim([0 1])
            xlim([0.5 nBin+0.5])
        end
    end
end

% Exp 5 - ScalingReprod NoMask

dx = [-0.5 0.5]*0.2;
c = 1;
yt = {-0.1:0.1:0.1,-0.2:0.2:0.2,[],[]};
tt = {'Scale','Reprod'};
for iseq = 1:2
    subplot(9,15,75+iseq)
    cla
    hold on
    for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
        iInd    = [1+(icond-1)*3 icond*3];
        x       = 1:2;
        icoefs  = squeeze(r{4}(iInd,iseq));
        ici     = squeeze(r_ci{4}(iInd,iseq,:));
        ipval   = squeeze(pval_fdr{4}(iInd,iseq));
        for ixrange = 1:length(x)
            plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:),'linewidth',lw)
        end
        plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
        ip = ipval<0.05;
        plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
    end
    plot([0.5 2.5],[0 0],'-k')
    xlim([0.5 2.5])
    if iseq == 1
        ylim([-0.2 0.2])
    else
        ylim([-0.3 0.3])
    end
    title(tt{iseq})
    set(gca,'fontsize',fz,'xtick',1:2,'ytick',yt{c},'xticklabel',[])
    c = c + 1;
end

% psychometric curve

nBin    = 4;
xval    = linspace(0,1,nBin+1);
for iseq = 1:2
    for ivar = 1:2
        for ig = 1:2
            subplot(18,16,190+2*iseq+16*(ivar-1)+ig)
            cla;
            hold on
            iInd    = results{4}.xy{iseq}.g==ig;
            iy      = results{4}.xy{iseq}.y(iInd);
            iyhat   = results{4}.xy{iseq}.yhat(iInd);
            switch ivar
                case 1
                    ix = results{4}.xy{iseq}.s0(iInd);
                case 2
                    ix = results{4}.xy{iseq}.s1(iInd);
            end
            ic1     = results{4}.xy{iseq}.d1(iInd);

            [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
            [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
            ip0         = order/size(ix,1);
            jx = NaN(size(ix));
            for iBin = 1:nBin
                if iBin == nBin
                    iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                else
                    iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                end
                jx(iInd) = iBin;
            end

            ing = length(unique(ic1));
            cl = jet(ing);
            for jc1 = 1:ing
                ky      = NaN(nBin,1);
                kyhat   = NaN(nBin,1);
                for iBin = 1:nBin
                    iInd        = (jx==iBin) & (ic1==jc1);
                    jy          = iy(iInd);
                    if ~isempty(jy)
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                end
                ix = find(~isnan(ky));
                plot(ix,kyhat(ix),'-','color',cl(jc1,:),'linewidth',1.2)
                plot(ix,ky(ix),'o','color','k','markerfacecolor',cl(jc1,:))
                if iseq == 1
                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                else
                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[-1 0 1],'fontsize',fz)
                end
            end
            if iseq == 1
                ylim([0 1])
            else
                ylim([-2 2])
            end
            xlim([0.5 nBin+0.5])
        end
    end
end

% Exp 6 - ScalingReprod

dx = [-0.5 0.5]*0.2;
c = 1;
yt = {-0.1:0.1:0.1,-0.2:0.2:0.2,[],[]};
tt = {'Scale','Reprod'};
for iseq = 1:2
    subplot(9,15,79+iseq)
    cla
    hold on
    for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
        iInd    = [1+(icond-1)*3 icond*3];
        x       = 1:2;
        icoefs  = squeeze(r{5}(iInd,iseq));
        ici     = squeeze(r_ci{5}(iInd,iseq,:));
        ipval   = squeeze(pval_fdr{5}(iInd,iseq));
        for ixrange = 1:length(x)
            plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:),'linewidth',lw)
        end
        plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
        ip = ipval<0.05;
        plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
    end
    plot([0.5 2.5],[0 0],'-k')
    xlim([0.5 2.5])
    if iseq == 1
        ylim([-0.2 0.2])
    else
        ylim([-0.3 0.3])
    end
    title(tt{iseq})
    set(gca,'fontsize',fz,'xtick',1:2,'ytick',yt{c},'xticklabel',[])
    c = c + 1;
end

% psychometric curve

nBin    = 4;
xval    = linspace(0,1,nBin+1);
for iseq = 1:2
    for ivar = 1:2
        for ig = 1:2
            subplot(18,16,194+2*iseq+16*(ivar-1)+ig)
            cla;
            hold on
            iInd    = results{5}.xy{iseq}.g==ig;
            iy      = results{5}.xy{iseq}.y(iInd);
            iyhat   = results{5}.xy{iseq}.yhat(iInd);
            switch ivar
                case 1
                    ix = results{5}.xy{iseq}.s0(iInd);
                case 2
                    ix = results{5}.xy{iseq}.s1(iInd);
            end
            ic1     = results{5}.xy{iseq}.d1(iInd);

            [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
            [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
            ip0         = order/size(ix,1);
            jx = NaN(size(ix));
            for iBin = 1:nBin
                if iBin == nBin
                    iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                else
                    iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                end
                jx(iInd) = iBin;
            end

            ing = length(unique(ic1));
            cl = jet(ing);
            for jc1 = 1:ing
                ky      = NaN(nBin,1);
                kyhat   = NaN(nBin,1);
                for iBin = 1:nBin
                    iInd        = (jx==iBin) & (ic1==jc1);
                    jy          = iy(iInd);
                    if ~isempty(jy)
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                end
                ix = find(~isnan(ky));
                plot(ix,kyhat(ix),'-','color',cl(jc1,:),'linewidth',1.2)
                plot(ix,ky(ix),'o','color','k','markerfacecolor',cl(jc1,:))
                if iseq == 1
                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                else
                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[-1 0 1],'fontsize',fz)
                end
            end
            if iseq == 1
                ylim([0 1])
            else
                ylim([-2 2])
            end
            xlim([0.5 nBin+0.5])
        end
    end
end

% Exp 7 - Motor

dx = [-0.5 0.5]*0.2;
yt = {-0.4:0.4:0.4,[],[],[]};
tt = {'Ring Same','Ring Diff','Pitch Same','Pitch Diff'};
for iseq = 1:4
    subplot(9,15,83+iseq)
    hold on
    for icond = 2:-1:1 % RT was not presented (to show RT; 4:-1:1)
        iInd    = [(1+(icond-1)*3) (icond*3)];
        x       = [1 2];
        icoefs  = squeeze(r{6}(iInd,iseq));
        ici     = squeeze(r_ci{6}(iInd,iseq,:));
        ipval   = squeeze(pval_fdr{6}(iInd,iseq));
        for ixrange = 1:length(x)
            plot([1 1]*x(ixrange)+dx(icond),ici(ixrange,:),'_-','color',colors(icond,:),'linewidth',lw)
        end
        plot(x+dx(icond),icoefs,'o','color',colors(icond,:),'MarkerFaceColor','w','linewidth',lw)
        ip = ipval<0.05;
        plot(x(ip)+dx(icond),icoefs(ip),'o','color',colors(icond,:),'MarkerFaceColor',colors(icond,:),'linewidth',lw)
    end
    plot([0.5 2.5],[0 0],'-k')
    xlim([0.5 2.5])
    ylim([-0.8 0.6])
    title(tt{iseq})
    set(gca,'fontsize',fz,'xtick',[1 2],'ytick',yt{iseq},'xticklabel',[])
end

% psychometric curve

nBin    = 3;
xval    = linspace(0,1,nBin+1);
ngs     = [2 4];
iInds   = [1 3 4 2];
c = 1;
for ifeature = 1:2
    for ihand = 1:2
        for ivar = 1:2
            for ig = 1:2
                subplot(18,16,200+2*(c-1)+16*(ivar-1)+ig)
                cla;
                hold on
                iInd    = results{6}.xy{ifeature,ihand}.g==ig;
                iy      = results{6}.xy{ifeature,ihand}.y(iInd);
                iyhat   = results{6}.xy{ifeature,ihand}.yhat(iInd);
                switch ivar
                    case 1
                        ix = results{6}.xy{ifeature,ihand}.s0(iInd);
                    case 2
                        ix = results{6}.xy{ifeature,ihand}.s1(iInd);
                end
                ic1     = results{6}.xy{ifeature,ihand}.d1(iInd);

                [~, rank]   = sort(ix, 1, 'ascend'); % 오름차순 정렬
                [~, order]  = sort(rank, 1, 'ascend'); % 정렬 순서 복원하여 순위 계산
                ip0         = order/size(ix,1);
                jx = NaN(size(ix));
                for iBin = 1:nBin
                    if iBin == nBin
                        iInd = (ip0>=xval(iBin)) & (ip0<=xval(iBin+1));
                    else
                        iInd = (ip0>=xval(iBin)) & (ip0<xval(iBin+1));
                    end
                    jx(iInd) = iBin;
                end

                ing = ngs(ig);
                cl = jet(ing);
                for jc1 = 1:ing
                    ky      = NaN(nBin,1);
                    kyhat   = NaN(nBin,1);
                    for iBin = 1:nBin
                        iInd        = (jx==iBin) & (ic1==jc1);
                        jy          = iy(iInd);
                        ky(iBin)    = mean(jy);
                        kyhat(iBin) = mean(iyhat(iInd));
                    end
                    plot(kyhat,'-','color',cl(jc1,:),'linewidth',1.2)
                    plot(ky,'o','color','k','markerfacecolor',cl(jc1,:))

                    set(gca,'xtick',1:nBin,'xticklabel',1:nBin,'ytick',[0 0.5 1],'fontsize',fz)
                end
                ylim([0 1])
                xlim([0.5 nBin+0.5])
            end
        end
        c = c + 1;
    end
end

%% Exp - ITI

tt = {'ITI short','ITI long'};
c = 1;
dx = linspace(-1,1,6)*0.2;
yt = {-0.7:0.7:0.7,[],[],[]};

ipExc = results{7}.results(end).pExcs;
iPExp = 0.14;
iInd0 = find(round(ipExc,2) == iPExp);
pvals = [];
betas = [];
for it = 1:2
    subplot(9,15,105+it) % stat1: t1s==1 & g1s==1; stat2: t1s==1 & g1s==2; stat3: t1s==2 & g1s==1; stat4: t1s==2 & g1s==2
    cla
    hold on
    for ic1s1 = 1:2
        if it == 1
            ibeta   = [results{7}.results(iInd0).stats{1}.beta(2*ic1s1); results{7}.results(iInd0).stats{2}.beta(2*ic1s1)]; % beta = [bias, d1, s0, s1, ...]
            ie      = 1.96*[results{7}.results(iInd0).stats{1}.se(2*ic1s1); results{7}.results(iInd0).stats{2}.se(2*ic1s1)];
            ih      = [fdr_bh(results{7}.results(iInd0).stats{1}.p) fdr_bh(results{7}.results(iInd0).stats{2}.p)];
            [~,~,~,ip] = fdr_bh(results{7}.results(iInd0).stats{1}.p);
            [~,~,~,jp] = fdr_bh(results{7}.results(iInd0).stats{2}.p);
            jbeta   = [results{7}.results(iInd0).stats{1}.beta results{7}.results(iInd0).stats{2}.beta];
        else
            ibeta   = [results{7}.results(iInd0).stats{3}.beta(2*ic1s1); results{7}.results(iInd0).stats{4}.beta(2*ic1s1)];
            ie      = 1.96*[results{7}.results(iInd0).stats{3}.se(2*ic1s1); results{7}.results(iInd0).stats{4}.se(2*ic1s1)];
            ih      = [fdr_bh(results{7}.results(iInd0).stats{3}.p) fdr_bh(results{7}.results(iInd0).stats{4}.p)];
            [~,~,~,ip] = fdr_bh(results{7}.results(iInd0).stats{3}.p);
            [~,~,~,jp] = fdr_bh(results{7}.results(iInd0).stats{4}.p);
            jbeta   = [results{7}.results(iInd0).stats{3}.beta results{7}.results(iInd0).stats{4}.beta];
        end
        if ic1s1 == 1
            betas   = cat(2,betas,jbeta);
            pvals   = cat(2,pvals,[ip jp]);
        end
        ih      = ih(2*ic1s1,:)';
        ici     = cat(2,ibeta-ie,ibeta+ie);
        
        ix = [1.25 1.75] + 3*(ic1s1-1);
        % plot(ix,ibeta,'o','color',colors(1+2*(it-1),:),'linewidth',lw)
        for i = 1:2
            plot([1 1]*ix(i),ici(i,:),'_-','color',colors(1+2*(i-1),:),'linewidth',lw)
            if ih(i) 
                plot(ix(i),ibeta(i),'o','color',colors(1+2*(i-1),:),'markerfacecolor',colors(1+2*(i-1),:),'linewidth',lw)
            else
                plot(ix(i),ibeta(i),'o','color',colors(1+2*(i-1),:),'markerfacecolor','w','linewidth',lw)
            end
        end
    end
    plot([0 6],[0 0],'k--')
    ylim([-1.2 1.2])
    xlim([0 6])
    title(tt{it})
    set(gca,'fontsize',fz,'xtick',[1.5 4.5],'ytick',yt{c},'xticklabel',[])
end

% Attribution probability

ng          = [2 4 8];
features    = [1 2];
nSub        = 58;
pars        = results{10}.pars;
for ifeature = features
    for ig = [1 2 3]
        ing = ng(ig);
        nStimBin = 8;
        Binval      = linspace(0,1,nStimBin+1);
        inTrans = NaN(ing,nStimBin,ing,nStimBin,nSub);
        for iSub = 1:nSub
            par         = pars{iSub,ifeature};
            iIndR       = find((par.condition == ig) & (par.StairTrainTest==3));
            d0          = par.Chc(1:45,iIndR);
            s0          = par.Stm(1:45,iIndR);
            inR         = length(iIndR);

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
            
            d1 = [NaN(1,inR); d0(1:end-1,:)];
            s1 = [NaN(1,inR); s0(1:end-1,:)];

            for id1 = 1:ing
                for is1 = 1:nStimBin
                    for is0 = 1:nStimBin
                        for id0 = 1:ing
                            iInd = (d0==id0) & (d1==id1) & (s1==is1) & (s0==is0);
                            inTrans(id0,is0,id1,is1,iSub) = sum(iInd(:));
                        end
                    end
                end
            end
        end

        inTrans0 = sum(inTrans,5);
        ifroms0 = [];
        ifromd1 = [];
        ifroms1 = [];
        for id0 = 1:ing
            ipTrans = squeeze(inTrans0(id0,:,:,:)); % s0,d1,s1

            ifroms0 = cat(1,ifroms0,mean(ipTrans./sum(ipTrans,1),[2 3],'omitnan')');
            ifromd1 = cat(1,ifromd1,mean(ipTrans./sum(ipTrans,2),[1 3],'omitnan'));
            ifroms1 = cat(1,ifroms1,squeeze(mean(ipTrans./sum(ipTrans,3),[1 2],'omitnan'))');
        end
        KL_s0 = sum(ifroms0 .* log(ifroms0 ./ (ones(size(ifroms0))./size(ifroms0,2))),2);
        KL_d10 = (sum(ifromd1 .* log(ifromd1 ./ (ones(size(ifromd1))./size(ifromd1,2))),2)./KL_s0)';
        KL_s10 = (sum(ifroms1 .* log(ifroms1 ./ (ones(size(ifroms1))./size(ifroms1,2))),2)./KL_s0)';
        
        subplot(18,16,248+1+16*(ig-1)+4*(ifeature-1))
        cla
        imagesc(ifroms0')
        clim([0 0.25])
        set(gca,'ytick',[1 nStimBin],'xtick',[1 ing],'fontsize',fz)

        subplot(18,16,248+2+16*(ig-1)+4*(ifeature-1))
        imagesc(ifromd1')
        set(gca,'ytick',[1 ing],'xtick',[1 ing],'fontsize',fz)
        if ig == 1
            clim([max([0 1/2-0.15]) 1/2+0.15])
        elseif ig == 2
            clim([max([0 1/4-0.15]) 1/4+0.15])
        elseif ig == 3
            clim([0 0.25])
        end
        
        subplot(18,16,248+3+16*(ig-1)+4*(ifeature-1))
        cla
        imagesc(ifroms1')
        clim([0 0.25])
        set(gca,'ytick',[1 nStimBin],'xtick',[1 ing],'fontsize',fz)

        subplot(18,16,248+4+16*(ig-1)+4*(ifeature-1))
        hold on
        plot(KL_d10,'-r','linewidth',1.2)
        plot(KL_s10,'-b','linewidth',1.2)
        plot([0 9],[1 1],'k--')
        ylim([0 2])
        xlim([0.5 0.5+length(KL_s10)])
        set(gca,'xtick',[1 ing],'ytick',[0 1 2],'fontsize',fz)
        
        if ig*ifeature == 1
            subplot(18,16,248+1-3)
            cb = colorbar;
            clim([0 0.25])
            cb.Ticks = [0, .1, 0.2]; % 원하는 tick 위치 설정
            set(gca,'fontsize',fz)
            subplot(18,16,248+1-2)
            cb = colorbar;
            clim([max([0 1/3-0.15]) 1/3+0.15])
            cb.Ticks = [0, .2, 0.4 0.6]; % 원하는 tick 위치 설정
            set(gca,'fontsize',fz)
            subplot(18,16,248+1-1)
            cb = colorbar;
            clim([max([0 1/2-0.15]) 1/2+0.15])
            cb.Ticks = [0, .2, 0.4 0.6]; % 원하는 tick 위치 설정
            set(gca,'fontsize',fz)
        end
    end
end
colormap(custom_colormap);

%

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5000 6000]/150)
saveas(gcf,'summary.png')