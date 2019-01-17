clear
close all
bbflagstring = ['N'; 'B'];
range_string = [4572; 9144; 18288];
phi_string = [0; 15; 30; 45];


for bbflagind = 1:length(bbflagstring)
for rangeind = 1:length(range_string)
for phiind = 1:length(phi_string)
filename = [bbflagstring(bbflagind) 'Bresults_range' num2str(range_string(rangeind)) '_phi' num2str(phi_string(phiind))];
load(filename);

%%
MCTind = 10;
if ~bbflag
%     % Plot single sample PSD
%     plotHPSD = Houtput.PSD(:,:,end);
%     [~, HDoAhat] = max(abs(plotHPSD(:)));  % finding estimated max DoA
%     [H_phiHatind, H_thetaHatind] = ind2sub(size(plotHPSD),HDoAhat);
%     
%     figure('position',[1950 10 1600 800])
%     Htxt = ['$\leftarrow \phi = $' num2str((Houtput.phi_estimate(MCTind)))...
%         '$^{\circ}, \theta = $' num2str(theta_search(H_thetaHatind))...
%         '$^{\circ}$, range = ' num2str(round(Houtput.range_estimate(MCTind)*10)/10) 'm'];
%     imagesc(phi_search,theta_search,10*log10(abs(plotHPSD)).')
%     set(gca,'Ydir','normal')
%     title(['H Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
%     xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
%     text(phi_search(H_phiHatind),theta_search(H_thetaHatind),Htxt,'FontSize',20)
%     hcb = colorbar; xlabel(hcb, 'dB')
%     caxis([-40 0]);
%     
%     % Plot averaged over all samples PSD
%     plotHPSD = mean(Houtput.PSD,3);
%     [~, HDoAhat] = max(abs(plotHPSD(:)));  % finding estimated max DoA
%     [H_phiHatind, H_thetaHatind] = ind2sub(size(plotHPSD),HDoAhat);
%     
%     figure('position',[1950 10 1600 800])
%     Htxt = ['$\leftarrow \phi = $' num2str((Houtput.phi_estimate(MCTind)))...
%         '$^{\circ}, \theta = $' num2str(theta_search(H_thetaHatind))...
%         '$^{\circ}$, range = ' num2str(round(Houtput.range_estimate(MCTind)*10)/10) 'm'];
%     imagesc(phi_search,theta_search,10*log10(abs(plotHPSD)).')
%     set(gca,'Ydir','normal')
%     title(['H$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
%     xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
%     text(phi_search(H_phiHatind),theta_search(H_thetaHatind),Htxt,'FontSize',20)
%     hcb = colorbar; xlabel(hcb, 'dB')
%     caxis([-40 0]);
%     
%     
%     plotVSPSD = VSoutput.PSD(:,:,end);
%     [~, VSDoAhat] = max(abs(plotVSPSD(:)));  % finding estimated max DoA
%     [VS_phiHatind, VS_thetaHatind] = ind2sub(size(plotVSPSD),VSDoAhat);
%     
%     figure('position',[1950 10 1600 800])
%     VStxt = ['$\leftarrow \phi = $' num2str((VSoutput.phi_estimate(MCTind)))...
%         '$^{\circ}, \theta = $' num2str(theta_search(VS_thetaHatind))...
%         '$^{\circ}$, range = ' num2str(round(VSoutput.range_estimate(MCTind)*10)/10) 'm'];
%     imagesc(phi_search,theta_search,10*log10(abs(mean(VSoutput.PSD,3))).')
%     set(gca,'Ydir','normal')
%     title(['VS Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) ' m'])
%     xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
%     text(phi_search(VS_phiHatind),theta_search(VS_thetaHatind),VStxt,'FontSize',20)
%     hcb = colorbar; xlabel(hcb, 'dB')
%     caxis([-40 0]);
%     
%     plotVSPSD = mean(VSoutput.PSD,3);
%     [~, VSDoAhat] = max(abs(plotVSPSD(:)));  % finding estimated max DoA
%     [VS_phiHatind, VS_thetaHatind] = ind2sub(size(plotVSPSD),VSDoAhat);
%     
%     figure('position',[1950 10 1600 800])
%     VStxt = ['$\leftarrow \phi = $' num2str((VSoutput.phi_estimate(MCTind)))...
%         '$^{\circ}, \theta = $' num2str(theta_search(VS_thetaHatind))...
%         '$^{\circ}$, range = ' num2str(round(VSoutput.range_estimate(MCTind)*10)/10) 'm'];
%     imagesc(phi_search,theta_search,10*log10(abs(mean(VSoutput.PSD,3))).')
%     set(gca,'Ydir','normal')
%     title(['VS$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) ' m'])
%     xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
%     text(phi_search(VS_phiHatind),theta_search(VS_thetaHatind),VStxt,'FontSize',20)
%     hcb = colorbar; xlabel(hcb, 'dB')
%     caxis([-40 0]);
    
    Houtput.mean_phi = mean(Houtput.phi_estimate);
    Houtput.std_phi = std(Houtput.phi_estimate);
    Houtput.mean_range = mean(Houtput.range_estimate);
    Houtput.std_range = std(Houtput.range_estimate);
    
    VSoutput.mean_phi = mean(VSoutput.phi_estimate);
    VSoutput.std_phi = std(VSoutput.phi_estimate);
    VSoutput.mean_range = mean(VSoutput.range_estimate);
    VSoutput.std_range = std(VSoutput.range_estimate);
    
    figure('position',[1950 10 1600 800])
    subplot(2,4,[1 2])
    hist(Houtput.phi_estimate,20);
    hold on
    %     plot([Houtput.mean_phi Houtput.mean_phi],[0 MCTs])
    %     errorbar(Houtput.mean_phi,MCTs/2,Houtput.std_phi)
    xlabel('$\hat{\phi}_H$ (deg)'); ylabel({'Occurences','$\times 10^4$'})
    title(['True DoA $\phi_s = $' num2str(ETS.az)...
        '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    
    
    subplot(2,4,[3 4])
    hist(Houtput.range_estimate)
    hold on
    %     plot([Houtput.mean_range Houtput.mean_range],[0 MCTs])
    %     errorbar(Houtput.mean_range,MCTs/2,Houtput.std_range)
    xlabel('$\hat{r}_H$ (m)'); ylabel({'$\times 10^4$'})
    title(['SNR$_{in}$ = ' num2str(SNRin) ' dB, '...
        num2str(MCTs) ' MCTs.'])
    
    
    subplot(2,4,[5 6])
    hist(VSoutput.phi_estimate)
    hold on
    %     plot([VSoutput.mean_phi VSoutput.mean_phi],[0 MCTs])
    %     errorbar(VSoutput.mean_phi,MCTs/2,VSoutput.std_phi)
    xlabel('$\hat{\phi}_{VS}$ (deg)'); ylabel({'Occurences','$\times 10^4$'})
    xlim([min(VSoutput.phi_estimate)-1 max(VSoutput.phi_estimate)+1])
    
    subplot(2,4,[7 8])
    hist(VSoutput.range_estimate)
    hold on
    %     plot([VSoutput.mean_range VSoutput.mean_range],[0 MCTs])
    %     errorbar(VSoutput.mean_range,MCTs/2,VSoutput.std_range)
    xlabel('$\hat{r}_{VS}$ (m)'); ylabel({'$\times 10^4$'})
    
else
    VSoutput.UIV.mean_phi = mean(VSoutput.UIV.phi_estimate,2);
    VSoutput.UIV.std_phi = std(VSoutput.UIV.phi_estimate,[],2);
    VSoutput.UIV.mean_range = mean(VSoutput.UIV.range_estimate);
    VSoutput.UIV.std_range = std(VSoutput.UIV.range_estimate);
    
    VSoutput.BB.mean_phi = mean(VSoutput.BB.phi_estimate,2);
    VSoutput.BB.std_phi = std(VSoutput.BB.phi_estimate,[],2);
    VSoutput.BB.mean_range = mean(VSoutput.BB.range_estimate);
    VSoutput.BB.std_range = std(VSoutput.BB.range_estimate);
    
    plotVSBBPSD = mean(VSoutput.BB.PSD,3);
    [~, BBDoAhat] = max(abs(plotVSBBPSD(:)));  % finding estimated max DoA
    [BBphiHatind, BBthetaHatind] = ind2sub(size(VSoutput.BB.PSD),BBDoAhat);
    
%     figure('position',[1950 10 1600 800])
%     BBtxt = ['$\leftarrow \phi = $' num2str(phi_search(BBphiHatind))...
%         '$^{\circ}, \theta = $' num2str(theta_search(BBthetaHatind))...
%         '$^{\circ}$, range = ' num2str(round(VSoutput.BB.range_estimate(MCTind)*10)/10) 'm'];
%     imagesc(phi_search,theta_search,10*log10(abs(plotVSBBPSD).'))
%     set(gca,'Ydir','normal')
%     title(['VS BB Beamforming, True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
%     xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
%     text(phi_search(BBphiHatind),theta_search(BBthetaHatind),BBtxt,'FontSize',20)
%     hcb = colorbar; xlabel(hcb, 'dB')
    
    figure('position',[1950 10 1600 800])
    subplot(2,4,[1 2])
    hold on
    for sensorind = 1:p.Nsensors
        hist(VSoutput.UIV.phi_estimate(sensorind,:))
    end
    %     plot([VSoutput.UIV.mean_phi VSoutput.UIV.mean_phi],[0 MCTs/p.Nsensors])
    %     errorbar(VSoutput.UIV.mean_phi,MCTs/2,VSoutput.UIV.std_phi)
    %     xlim([min(VSoutput.UIV.phi_estimate)-1 max(VSoutput.UIV.phi_estimate)+1])
    xlabel('$\hat{\phi}_{uiv}$'); ylabel({'Occurences','$\times 10^4$'})
    title(['SNR$_{in}$ = ' num2str(SNRin) ' dB,'...
        num2str(MCTs) ' MCTs.'])
    
    subplot(2,4,[3 4])
    hist(VSoutput.UIV.range_estimate)
    hold on
    xlim([min(VSoutput.UIV.range_estimate)-1 max(VSoutput.UIV.range_estimate)+1])
    %     plot([VSoutput.UIV.mean_range VSoutput.UIV.mean_range],[0 MCTs/p.Nsensors])
    %     errorbar(VSoutput.UIV.mean_range,MCTs/2,VSoutput.UIV.std_range)
    title([' True DoA: $\phi_s = $' num2str(ETS.az)...
        '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    xlabel('$\hat{r}_{uiv}$'); ylabel({'$\times 10^4$'})
    
    subplot(2,4,[5 6])
    hold on
    for sensorind = 1:p.Nsensors
        hist(VSoutput.BB.phi_estimate(sensorind,:))
    end
    %     plot([VSoutput.BB.mean_phi VSoutput.BB.mean_phi],[0 MCTs/p.Nsensors])
    %     errorbar(VSoutput.BB.mean_phi,MCTs/2,VSoutput.BB.std_phi)
    xlabel('$\hat{\phi}_{BB}$'); ylabel({'Occurences','$\times 10^4$'})
    xlim([min(mean(VSoutput.BB.phi_estimate,1))-1 max(mean(VSoutput.BB.phi_estimate,1))+1])
    
    subplot(2,4,[7 8])
    hist(VSoutput.BB.range_estimate)
    hold on
    %     plot([VSoutput.BB.range_estimate VSoutput.BB.range_estimate],[0 MCTs/p.Nsensors])
    %     xlim([min(VSoutput.UIV.range_estimate)-5 max(VSoutput.UIV.range_estimate)+5])
    %     errorbar(VSoutput.BB.mean_range,MCTs/2,VSoutput.BB.std_range)
    xlabel('$\hat{r}_{BB}$'); ylabel({'$\times 10^4$'})
    
end

end
end
end

