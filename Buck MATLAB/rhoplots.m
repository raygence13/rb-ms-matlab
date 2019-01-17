%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    rho vs alpha plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load  rhoanalysis.mat
f10 = figure(10);
set(f10,'DefaultAxesFontSize',24)
for b = 1:1
    clf
    semilogx(1-alpha,rhoNoiseder,'r','LineWidth',4)
    hold on
    semilogx(1-alpha(1:2:end),abs(rhoNoisesim(1:2:end)),':dr','MarkerSize',9)
%     semilogx(1-alpha,(rhoSigNoiseder(b,:)),'LineWidth',4,'color',[0 0 1])
%     semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim(1:2:end)),':d','MarkerSize',9,'color',[0 0 1])
    semilogx(1-alpha,(rhoSigNoiseder_20dB),'m','LineWidth',4)%,'color',[1 .6 0])
    semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_20dB(1:2:end)),'m:d','MarkerSize',9)%,'color',[1 .6 0])
    semilogx(1-alpha,(rhoSigNoiseder_10dB),'b','LineWidth',4)
    semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_10dB(1:2:end)),'b:d','MarkerSize',9)
    semilogx(1-alpha,(rhoSigNoiseder_0dB),'LineWidth',4,'color',[0 .6 .1])
    semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_0dB(1:2:end)),':d','MarkerSize',9,'color',[0 .6 .1])
    grid on
    ax = gca;
    xtick = linspace(.1,1,10);
    set(ax,'XTick',xtick);
    xticklabel = abs(round(xtick.*100)./100 - 1);
    set(ax,'XTickLabel',xticklabel);
    set(ax,'Xdir','reverse')
%     ax.XTickLabelRotation=45;
    xlim([0.01 1]); ylim([0 1.1]);
    xlabel('$\alpha$')
    ylabel('$\rho$')
    legend('$\rho_{0}$','','$\rho_{1}$','SNR$_{in} = -20$dB','$\rho_{1}$','SNR$_{in} = -10$dB','$\rho_{1}$','SNR$_{in} = 0$dB','Location','SouthEast')
    title(['Correlation Coefficient $\rho$ vs. $\alpha$, M = ' num2str(M) ', N = ' num2str(N)])
    
end