function [ VSoutput, Houtput ] = ProcessETS(bbflag, SNRin, phi_vec, theta_vec, MCTs, ETS, p)
%[ VSoutput, Houtput ] = ProcessETS(bbflag, ETS, p, SNRin, phi_vec, theta_vec, MonteCarloTrials)
% Function to process element time series
% RB 9MAR2018. RB 15MAR2018.
% Inputs:
% bbflag:       broadband signal option
% SNRin:        sensor level SNR
% phi_vec:      user defined phi vector (deg)
% theta_vec:    user defined theta vector (deg)
% MCTs:         number of Monte-Carlo Trials to perform
%
% ETS:          structure containing h,x,y,z ETS
% ETS.h:        hydrophone element time series [Nfft x Nsensors] [Pa](t)
% ETS.Hs:       hydrophone sensitivity [dB re V/uPa]
% ETS.a:        accelerometer element time series [m/s^2](t)
% ETS.x:        x-component element time series [Nfft x Nsensors] [m/s^2](t)
% ETS.Axs:      x sensitivity [V/g]
% ETS.y:        y-component element time series [Nfft x Nsensors] [m/s^2](t)
% ETS.Ays:      y sensitivity [V/g]
% ETS.z:        z-component element time series [Nfft x Nsensors] [m/s^2](t)
% ETS.Azs:      z sensitivity [V/g]
% ETS.time:     time vector [1 x Nfft] [s]
% ETS.fs:       sampling frequency [Hz]
% ETS.Nfft:     fft size
% ETS.f0:       center frequency [Hz]
% ETS.az:       target azimuth bearing [deg]
% ETS.de:       target d/e bearing [deg]
% ETS.range:    target range [m]
%
% p:                structure containing general parameters
% p.AEL:            Array Element Locations [3 x Nsensors] [m]
% p.Nsensors:       number of sensors
% p.MaxDistance:    max distance between sensor pairs in geometry [m]
% p.depth:          depth of array [m]
% p.c:              sound speed of water [m/s]
% p.rho:            density [kg/m^3]
% p.f0:             center frequency
%
% Outputs:
% VSoutput:     structure containing VS outputs
% Houtput:     structure containing H outputs

if max(abs(theta_vec)) <= 1
    theta_vec = asind(theta_vec);
end
l_phi = length(phi_vec);
l_theta = length(theta_vec);

% steering vectors
ux = cosd(phi_vec)'*cosd(theta_vec);
ux = ux(:);     % x-component of unit vector
uy = sind(phi_vec)'*cosd(theta_vec);
uy = uy(:);     % z-component of unit vector
uz = ones(l_phi,1)*sind(theta_vec);
uz = uz(:);     % y-component of unit vector

% 3-dimensional unit vector
u = zeros(3,l_phi*l_theta);
u(1,:) = ux;
u(2,:) = uy;
u(3,:) = uz;

varP = trace(ETS.h*ETS.h')/p.Nsensors;
varNp = varP/(10^(SNRin/10));
varA = trace(ETS.a*ETS.a')/p.Nsensors;
varNa = varA/(10^(SNRin/10));

% initialize memory
if ~bbflag
    Houtput.PSD = zeros(l_phi,l_theta,MCTs);
    Houtput.phi_estimate = zeros(1,MCTs);
    Houtput.range_estimate = zeros(1,MCTs);
    VSoutput.PSD = zeros(l_phi,l_theta,MCTs);
    VSoutput.phi_estimate = zeros(1,MCTs);
    VSoutput.range_estimate = zeros(1,MCTs);
else
    Houtput = [];
    VSoutput.UIV.phi_estimate = zeros(p.Nsensors,MCTs);
    VSoutput.UIV.range_estimate = zeros(1,MCTs);
    VSoutput.BB.PSD = zeros(l_phi,l_theta,MCTs);
    VSoutput.BB.phi_estimate = zeros(p.Nsensors,MCTs);
    VSoutput.BB.range_estimate = zeros(1,MCTs);
end

for MCTind = 1:MCTs
    
    % Add noise
    hData = ETS.h + varNp.*randn(ETS.Nfft,p.Nsensors);
    xData = ETS.x + varNa.*randn(ETS.Nfft,p.Nsensors);
    yData = ETS.y + varNa.*randn(ETS.Nfft,p.Nsensors);
    zData = ETS.z + varNa.*randn(ETS.Nfft,p.Nsensors);
    
    %% Beamforming
    HData = fft(hData/10^(ETS.Hs/20)*1e-6,ETS.Nfft,1);
    % hydrophone spectrum converted from v(t) to Pa(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    
    XData = fft(xData*9.8/ETS.Axs,ETS.Nfft,1);
    % X acceleration spectrum converted from v(t) to m/s^2(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    X_p = bsxfun(@times,XData*p.rho*p.c,1./(1i*2*pi*ETS.frequencies.'));
    % X pressure spectrum converted from m/s^2(f) to Pa(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    YData = fft(yData*9.8/ETS.Ays,ETS.Nfft,1);
    % Y acceleration spectrum converted from v(t) to m/s^2(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    Y_p = bsxfun(@times,YData*p.rho*p.c,1./(1i*2*pi*ETS.frequencies.'));
    % Y pressure spectrum converted from m/s^2(f) to Pa(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    ZData = fft(zData*9.8/ETS.Azs,ETS.Nfft,1);
    % Z acceleration spectrum converted from v(t) to m/s^2(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    Z_p = bsxfun(@times,ZData*p.rho*p.c,1./(1i*2*pi*ETS.frequencies.'));
    % Z pressure spectrum converted from m/s^2(f) to Pa(f) [ETS.Nfft x p.Nsensors x Nsnapshots]
    
    if ~bbflag
        f0idx = round(p.f0/ETS.fs*ETS.Nfft)+1;                    % bin index for f0
        % Hydrophone only processing
        H_SM    = exp(-1i*2*pi/p.c*p.f0*p.AEL.'*u)/p.Nsensors;   % Hydrophone steering matrix [p.Nsensors x (l_phi*l_theta)]
        H_snaps = HData(f0idx,:).';    % Hydrophone Nsnapshots [1 x p.Nsensors]
        H_BR    = mean(H_snaps'*H_SM,1); H_BR = reshape(H_BR,l_phi,l_theta);
        % Hydrophone beamformed response averaged over Nsnapshots and reshaped [l_phi x l_theta]
        H_PSD   = H_BR.*conj(H_BR);           % Power spectral density [l_phi x l_theta]
        Houtput.PSD(:,:,MCTind)   = H_PSD./max(max(H_PSD));     % Normalized PSD
        
        [~, HDoAhat] = max(abs(H_PSD(:)));  % finding estimated max DoA
        [H_phiHatind, H_thetaHatind] = ind2sub(size(H_PSD),HDoAhat);
        Houtput.phi_estimate(MCTind) = phi_vec(H_phiHatind); H_thetaHat = theta_vec(H_thetaHatind);
        Houtput.range_estimate(MCTind) = -p.depth/sind((H_thetaHat));
        Houtput.range_estimate(isinf(Houtput.range_estimate)) = nan;
        
        % VS steering vectors [p.Nsensors x (l_phi*l_theta)]
        X_SM = bsxfun(@times,H_SM,ux.');
        Y_SM = bsxfun(@times,H_SM,uy.');
        Z_SM = bsxfun(@times,H_SM,uz.');
        % Accelerometer steering matrices
        alpha = 0.5;
        % alpha=0.5 for cardioid
        VS_W = [H_SM*alpha;X_SM*(1-alpha);Y_SM*(1-alpha);Z_SM*(1-alpha)]/sqrt(p.Nsensors);
        % combined steering matrices
        
        % XYZ snapshots [p.Nsensors x Nsnapshots]
        X_snaps = X_p(f0idx,:).';
        Y_snaps = Y_p(f0idx,:).';
        Z_snaps = Z_p(f0idx,:).';
        
        VS_BR = mean([H_snaps;X_snaps;Y_snaps;Z_snaps]'*VS_W,1);
        % VS beamformed response averaged over Nsnapshots [1 x (l_phi*l_theta)]
        VS_BR = reshape(VS_BR,l_phi,l_theta);
        % Reshaped beamformed response [l_phi x l_theta]
        VS_PSD = VS_BR.*conj(VS_BR);
        VSoutput.PSD(:,:,MCTind) = VS_PSD./max(max(VS_PSD));
        
        [~, VSDoAhat] = max(abs(VS_PSD(:)));  % finding estimated max DoA
        [VS_phiHatind, VS_thetaHatind] = ind2sub(size(VS_PSD),VSDoAhat);
        VSoutput.phi_estimate(MCTind) = phi_vec(VS_phiHatind); VS_thetaHat = theta_vec(VS_thetaHatind);
        VSoutput.range_estimate(MCTind) = -p.depth/sind((VS_thetaHat));
        
        
    else
        Ha = zeros(size(HData));    % Initialize hydrophone acceleration spectrum
        % Convert hydrophone pressure data to acceleration
        Ha(1:ETS.Nfft/2,:,:) = bsxfun(@times,HData(1:ETS.Nfft/2,:,:)./(p.rho*p.c),(1i*2*pi.*ETS.frequencies(1:ETS.Nfft/2).'));
        Ha(ETS.Nfft/2+2:ETS.Nfft,:,:) = bsxfun(@times,HData(ETS.Nfft/2+2:ETS.Nfft,:,:)./(p.rho*p.c),(-1i*2*pi.*ETS.frequencies(2:ETS.Nfft/2).'));
        
        Ah = ifft(Ha,ETS.Nfft,1);   % hydrophone acceleration time series [ETS.Nfft x p.Nsensors]
        Ax = xData*9.8/ETS.Axs;
        Ay = yData*9.8/ETS.Ays;
        Az = zData*9.8/ETS.Azs;
        
        % Calculation of bearing using intensity at each sensor
        uiv = zeros(3,p.Nsensors);    % initialize memory for unit intensity vector
        for sensorInd = 1:p.Nsensors
            s = real([Ax(:,sensorInd) Ay(:,sensorInd) Az(:,sensorInd)]'*Ah(:,sensorInd));
            %   intensity vector calculated from product of Axyz with Ah
            uiv(:,sensorInd) = s./(sqrt(sum(s.*s)));
        end
        VSoutput.UIV.phi_estimate(:,MCTind) = atand(uiv(2,:)./uiv(1,:));
        VS_UIV_thetaHat = asind(uiv(3,:));
        VSoutput.UIV.range_estimate(MCTind) = abs(p.depth/sind(mean(VS_UIV_thetaHat)));
        
        % Calculation of bearing using conventional beamforming
        BB_BR = zeros(l_phi*l_theta,p.Nsensors);  % Initialize memory for broadband PSD
        for sensorInd = 1:p.Nsensors              % Do broadband beamforming at each sensor
            Ahxyz = [Ah(:,sensorInd) Ax(:,sensorInd) Ay(:,sensorInd) Az(:,sensorInd)]'*[Ah(:,sensorInd) Ax(:,sensorInd) Ay(:,sensorInd) Az(:,sensorInd)];
            % [4 x 4] correlation matrix
            for beamind = 1:(l_phi*l_theta)     % Loop over beams
                BB_BR(beamind,sensorInd) = [1;u(:,beamind)]'*Ahxyz*[1;u(:,beamind)];
            end
            BB_BR(:,sensorInd) = BB_BR(:,sensorInd)/max(abs(BB_BR(:,sensorInd)));
            % normalizing max value to 1
        end
        VSoutput.BB.PSD(:,:,MCTind) = reshape(mean(BB_BR,2),l_phi,l_theta);
        temp = VSoutput.BB.PSD(:,:,MCTind);
        [~, BBDoAhat] = max(abs(temp(:)));  % finding estimated max DoA
        [BBphiHatind, BBthetaHatind] = ind2sub(size(temp),BBDoAhat);
        VSoutput.BB.phi_estimate(:,MCTind) = phi_vec(BBphiHatind); VS_BB_thetaHat = theta_vec(BBthetaHatind);
        VSoutput.BB.range_estimate(MCTind) = abs(p.depth/sind(mean(VS_BB_thetaHat)));
        
    end
end


if ~bbflag
    % Plot single sample PSD
    plotHPSD = Houtput.PSD(:,:,end);
    [~, HDoAhat] = max(abs(plotHPSD(:)));  % finding estimated max DoA
    [H_phiHatind, H_thetaHatind] = ind2sub(size(plotHPSD),HDoAhat);
    
    figure('position',[1950 10 1600 800])
    Htxt = ['$\leftarrow \phi = $' num2str((Houtput.phi_estimate(MCTind)))...
        '$^{\circ}, \theta = $' num2str(theta_vec(H_thetaHatind))...
        '$^{\circ}$, range = ' num2str(round(Houtput.range_estimate(MCTind)*10)/10) 'm'];
    imagesc(phi_vec,theta_vec,10*log10(abs(plotHPSD)).')
    set(gca,'Ydir','normal')
    title(['H$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
    text(phi_vec(H_phiHatind),theta_vec(H_thetaHatind),Htxt,'FontSize',20)
    hcb = colorbar; xlabel(hcb, 'dB')
        caxis([-40 0]);
        
    % Plot averaged over all samples PSD
    plotHPSD = mean(Houtput.PSD,3);
    [~, HDoAhat] = max(abs(plotHPSD(:)));  % finding estimated max DoA
    [H_phiHatind, H_thetaHatind] = ind2sub(size(plotHPSD),HDoAhat);
    
    figure('position',[1950 10 1600 800])
    Htxt = ['$\leftarrow \phi = $' num2str((Houtput.phi_estimate(MCTind)))...
        '$^{\circ}, \theta = $' num2str(theta_vec(H_thetaHatind))...
        '$^{\circ}$, range = ' num2str(round(Houtput.range_estimate(MCTind)*10)/10) 'm'];
    imagesc(phi_vec,theta_vec,10*log10(abs(plotHPSD)).')
    set(gca,'Ydir','normal')
    title(['H$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
    text(phi_vec(H_phiHatind),theta_vec(H_thetaHatind),Htxt,'FontSize',20)
    hcb = colorbar; xlabel(hcb, 'dB')
        caxis([-40 0]);
    
    
    plotVSPSD = VSoutput.PSD(:,:,end);
    [~, VSDoAhat] = max(abs(plotVSPSD(:)));  % finding estimated max DoA
    [VS_phiHatind, VS_thetaHatind] = ind2sub(size(plotVSPSD),VSDoAhat);
    
    figure('position',[1950 10 1600 800])
    VStxt = ['$\leftarrow \phi = $' num2str((VSoutput.phi_estimate(MCTind)))...
        '$^{\circ}, \theta = $' num2str((VS_thetaHat))...
        '$^{\circ}$, range = ' num2str(round(VSoutput.range_estimate(MCTind)*10)/10) 'm'];
    imagesc(phi_vec,theta_vec,10*log10(abs(mean(VSoutput.PSD,3))).')
    set(gca,'Ydir','normal')
    title(['VS$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) ' m'])
    xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
    text(phi_vec(VS_phiHatind),theta_vec(VS_thetaHatind),VStxt,'FontSize',20)
    hcb = colorbar; xlabel(hcb, 'dB')
        caxis([-40 0]);
        
    plotVSPSD = mean(VSoutput.PSD,3);
    [~, VSDoAhat] = max(abs(plotVSPSD(:)));  % finding estimated max DoA
    [VS_phiHatind, VS_thetaHatind] = ind2sub(size(plotVSPSD),VSDoAhat);
    
    figure('position',[1950 10 1600 800])
    VStxt = ['$\leftarrow \phi = $' num2str((VSoutput.phi_estimate(MCTind)))...
        '$^{\circ}, \theta = $' num2str((VS_thetaHat))...
        '$^{\circ}$, range = ' num2str(round(VSoutput.range_estimate(MCTind)*10)/10) 'm'];
    imagesc(phi_vec,theta_vec,10*log10(abs(mean(VSoutput.PSD,3))).')
    set(gca,'Ydir','normal')
    title(['VS$_{avg}$ Response, SNR$_{in}$ = ' num2str(SNRin) ' dB. True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) ' m'])
    xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
    text(phi_vec(VS_phiHatind),theta_vec(VS_thetaHatind),VStxt,'FontSize',20)
    hcb = colorbar; xlabel(hcb, 'dB')
        caxis([-40 0]);
    
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
    xlabel('$\hat{\phi}_H$')
    title(['SNRin$_{in}$ = ' num2str(SNRin) ' dB, '...
        num2str(MCTs) ' MCTs.'])
    
    
    subplot(2,4,[3 4])
    hist(Houtput.range_estimate)
    hold on
%     plot([Houtput.mean_range Houtput.mean_range],[0 MCTs])
%     errorbar(Houtput.mean_range,MCTs/2,Houtput.std_range)
    xlabel('$\hat{r}_H$')
    title(['True DoA $\phi_s = $' num2str(ETS.az)...
        '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    
    subplot(2,4,[5 6])
    hist(VSoutput.phi_estimate)
    hold on
%     plot([VSoutput.mean_phi VSoutput.mean_phi],[0 MCTs])
%     errorbar(VSoutput.mean_phi,MCTs/2,VSoutput.std_phi)
    xlabel('$\hat{\phi}_{VS}$')
    xlim([min(VSoutput.phi_estimate)-1 max(VSoutput.phi_estimate)+1])
    
    subplot(2,4,[7 8])
    hist(VSoutput.range_estimate)
    hold on
%     plot([VSoutput.mean_range VSoutput.mean_range],[0 MCTs])
%     errorbar(VSoutput.mean_range,MCTs/2,VSoutput.std_range)
    xlabel('$\hat{r}_{VS}$')
    
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
    
    figure('position',[1950 10 1600 800])
    BBtxt = ['$\leftarrow \phi = $' num2str(phi_vec(BBphiHatind))...
        '$^{\circ}, \theta = $' num2str(theta_vec(BBthetaHatind))...
        '$^{\circ}$, range = ' num2str(round(VSoutput.BB.range_estimate(MCTind)*10)/10) 'm'];
    imagesc(phi_vec,theta_vec,10*log10(abs(plotVSBBPSD).'))
    set(gca,'Ydir','normal')
    title(['VS BB Beamforming, True DoA $\phi_s = $' num2str(ETS.az) '$^{\circ}, \theta_s = $' num2str(round(ETS.de*10)/10) '$^{\circ}$ and range = ' num2str(ETS.range) 'm'])
    xlabel('Port $\leftarrow \phi \rightarrow$ Stbd'); ylabel('$\theta$, in deg $\rightarrow$ Up')
    text(phi_vec(BBphiHatind),theta_vec(BBthetaHatind),BBtxt,'FontSize',20)
    hcb = colorbar; xlabel(hcb, 'dB')
    
    figure('position',[1950 10 1600 800])
    subplot(2,4,[1 2])
    hold on
    for sensorind = 1:p.Nsensors
        hist(VSoutput.UIV.phi_estimate(sensorind,:))
    end
%     plot([VSoutput.UIV.mean_phi VSoutput.UIV.mean_phi],[0 MCTs/p.Nsensors])
%     errorbar(VSoutput.UIV.mean_phi,MCTs/2,VSoutput.UIV.std_phi)
%     xlim([min(VSoutput.UIV.phi_estimate)-1 max(VSoutput.UIV.phi_estimate)+1])
    xlabel('$\hat{\phi}_{uiv}$')
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
    xlabel('$\hat{r}_{uiv}$')
    
    subplot(2,4,[5 6])
    hold on
    for sensorind = 1:p.Nsensors
        hist(VSoutput.BB.phi_estimate(sensorind,:))
    end
%     plot([VSoutput.BB.mean_phi VSoutput.BB.mean_phi],[0 MCTs/p.Nsensors])
%     errorbar(VSoutput.BB.mean_phi,MCTs/2,VSoutput.BB.std_phi)
    xlabel('$\hat{\phi}_{BB}$')
    xlim([min(mean(VSoutput.BB.phi_estimate,1))-1 max(mean(VSoutput.BB.phi_estimate,1))+1])
    
    subplot(2,4,[7 8])
    hist(VSoutput.BB.range_estimate)
    hold on
%     plot([VSoutput.BB.range_estimate VSoutput.BB.range_estimate],[0 MCTs/p.Nsensors])
%     xlim([min(VSoutput.UIV.range_estimate)-5 max(VSoutput.UIV.range_estimate)+5])
%     errorbar(VSoutput.BB.mean_range,MCTs/2,VSoutput.BB.std_range)
    xlabel('$\hat{r}_{BB}$')
    
end
end

