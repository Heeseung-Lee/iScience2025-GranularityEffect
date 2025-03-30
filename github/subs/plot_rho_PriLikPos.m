function plot_rho_PriLikPos(isplot,iT,nT,rangetheta,itheta_lik,itheta_pos,itheta_pri0,itheta_pri,itheta_pos_conded,...
    rangerho0,rangerho,irho_lik0,irho_lik,irho_lik00,irho_pos,irho_pri0,irho_pri,irho_pos_conded)

if isplot == 1 || (isplot == 2 && iT == nT)

    figure(1)
    clf
    subplot(2,1,1)
    hold on
    plot(rangetheta,itheta_lik,'k-','linewidth',1.5)
    plot(rangetheta,itheta_pri0,'b--','linewidth',1.5)
    plot(rangetheta,itheta_pos,'r--','linewidth',1.5)
    plot(rangetheta,itheta_pos_conded,'r-','linewidth',1.5)
    plot(rangetheta,itheta_pri,'b-','linewidth',1.5)
    xlim([-3 3])
    xlabel('theta')
    ylabel('probability')
    set(gca,'fontsize',15)
    legend({'lik','prior-pre','pos','pos-conditioned','pri-updated'},'Location','eastoutside')

    subplot(2,1,2)
    hold on
    plot(rangerho0,irho_lik0,'--','color',[1 1 1]*0.6,'linewidth',1.5)
    plot(rangerho,irho_lik,'k-','linewidth',1.5)
    plot(rangerho,irho_lik-irho_lik00,'k:','color',[1 1 1]*0.6,'linewidth',1.5)
    plot(rangerho,irho_pri0,'b--','linewidth',1.5)
    plot(rangerho,irho_pos,'r--','linewidth',1.5)
    plot(rangerho,irho_pos_conded,'r-','linewidth',1.5)
    plot(rangerho,irho_pri,'b-','linewidth',1.5)
    
    grid on
    xlim([-0.2 1.2])
    xlabel('rho')
    ylabel('probability')
    set(gca,'fontsize',15)
    legend({'lik-pre','lik-concat','lik-concat - lik-pre','prior-pre',...
        'pos','pos-condtioned','prior-updated'},'Location','eastoutside')
    set(gcf, 'Position', [1000 1500 1200 900])

end

end