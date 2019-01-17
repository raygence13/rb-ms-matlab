close all
clear
present(0)
%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debugflag = 0;                  % iterates only one block of data if 1
do_VSP = 1;

headersize  = 104;
Fs          = 12500;        	% Sampling Frequency
Nchans      = 384;             	% H + X + Y + Z
Nsens       = Nchans/4;         % number of H and VS sensors
c           = 1500;            	% medium sound speed [m/s]
rho         = 1026;          	% medium density [kg/m^3]

d           = 0.1406525;        % channel spacing [m]
idx_h       = 1:96;             % 48 fwd VS elements, 48 aft VS elements
idx_x       = idx_h + 96;
idx_y       = idx_x + 96;
idx_z       = idx_y + 96;

Aperture    = d*96*10;           % Overestimated aperture length [m]
prop_time   = round(Aperture/c*Fs*10);       % propagation time [bins]
nPts = prop_time;               % block size [bins]

p = [zeros(Nsens,1)...       % sensors along z-axis
    zeros(Nsens,1)...
    (0:Nsens-1)'*d].';       % sensor position vector [x y z] 3 x Nsens

load('FWD_Calibration.mat');
% f0  = 1000:25:6000;          % frequency of interest [Hz]
f0 = FWD_Cal_H(1,1:end).*1e3;
FWD_Cal_H(2:end,:) = FWD_Cal_H(2:end,:)./1e6;
l_f0 = length(f0);

tfft  = 2^nextpow2(2*Fs);       % time series fft size
bins  = 1:tfft/2;               % bin number for first half of frequency vector
f0idx = floor(f0.*tfft./Fs);    % frequency of interest index [bin#]

time  = ((1:nPts)*(1/Fs));      % time vector for block data
f     = Fs*(0:tfft/2-1)/tfft;   % frequency vector for fft
overlap = 1/4;                   % block overlap percentage
OverLapPts = round(nPts*(1-overlap));   % overlap size
h_window = hanning(nPts);       % time series window
blocks_to_avg = round(1/overlap);

% 3-dimensional beamforming parameters
phi = -180:180;      % degrees
l_phi = length(phi);
theta = 88:2:100;
l_theta = length(theta);

ux = repmat(cosd(theta),[1 l_phi]);   % sensors along z-axis
uy = sind(phi)'*sind(theta);
uz = cosd(phi)'*sind(theta);

ux = ux(:);     % x-component of unit vector
uy = uy(:);     % z-component of unit vector
uz = uz(:);     % y-component of unit vector

% 3-dimensional unit vector
u = zeros(3,l_phi*l_theta);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);

% 1-dimensional beamforming parameter
ucbf = cosd(phi);

% [fname, pname] = uigetfile('*.rst', 'select input RST data file');
[fname, pname] = uigetfile('*.ets', 'select input ETS data file');
if(~ischar(fname))
    return;
end
[nullpath, filename, theext]=fileparts(fname);
cd(pname);
input_fullname = [pname fname];
fidin    = fopen(input_fullname,'r','l');
[VSData,RollData,SamplesToRead] = GetFLVSTASamples( fidin, nPts );

if debugflag
    numBlocks = 50;
else
    numBlocks = floor(SamplesToRead/nPts);
%     numBlocks = 100;
end

H_DEvsAZvsF = zeros(l_phi,l_theta,l_f0,numBlocks);
VS_DEvsAZvsF = zeros(l_phi,l_theta,l_f0,numBlocks);

profile on
%%
for block = 1:numBlocks
    
    VSData = VSData - repmat(mean(VSData,2),1,nPts); % Remove Mean
    blockdata = bsxfun(@times,VSData,h_window.');
    
    % Parsing out hydrophone and accelerometer time series
    % Computing temporal frequencies for each channel
    h   = blockdata(idx_h,:);
    % hydrophone time series
    H  = fft(h,tfft,2)./tfft; H = H(:,bins);
    % H pressure spectrum
    
    if do_VSP
    x = blockdata(idx_x,:);
    % x time series
    X = fft(x,tfft,2)./tfft; X = X(:,bins);
    % X acceleration spectrum
    X_p = bsxfun(@times,X*rho*c,1./(1i*2*pi*f));
    % X pressure spectrum
    
    y = blockdata(idx_y,:);
    % y time series
    Y = fft(y,tfft,2)./tfft; Y = Y(:,bins);
    % Y acceleration spectrum
    Y_p = bsxfun(@times,Y*rho*c,1./(1i*2*pi*f));
    % Y pressure spectrum
    
    z = blockdata(idx_z,:);
    % z time series
    Z = fft(z,tfft,2)./tfft; Z = Z(:,bins);
    % Z acceleration spectrum
    Z_p = bsxfun(@times,Z*rho*c,1./(1i*2*pi*f));
    % Z pressure spectrum
    
    RollAvg = mean(RollData,2);
    Rotation_z = sind(RollAvg);
    Rotation_y = zeros(96,1);
    Rotation_x = cosd(RollAvg);
    % Avg over time roll measurements
    end
    
    % beamforming data
    for indf = 1:l_f0 % Outputs are H/VS_DEvsAZ, H/VS_FRAZ
        
        % Hydrophone only processing
        H_SM    = exp(-1i*2*pi/c*f0(indf)*p.'*u);
        % Hydrophone steering matrix
        H_snaps = H(:,f0idx(indf)); Hn_snaps = H_snaps/max(abs(H_snaps));
        % Hydrophone snapshots
        H_BR = H_SM'*H_snaps/96; H_BR = reshape(H_BR,l_phi,l_theta);
        % Hydrophone beam response
        H_DEvsAZvsF(:,:,indf,block) = H_BR.*conj(H_BR);
        % Hyrdophone DE vs. AZ PSD
        
        if do_VSP
        % VS processing
        X_SM = H_SM.*(Rotation_x*ux.');
        Y_SM = H_SM.*(Rotation_z*uy.');
        Z_SM = H_SM.*(Rotation_z*uz.');
        % Accelerometer steering matrices
        alpha = 0.5;
        % alpha=0.5 for cardioud
        VS_W = [H_SM*alpha;X_SM*(1-alpha);Y_SM*(1-alpha);Z_SM*(1-alpha)];
        % combined steering matrices
        
        X_snaps = X_p(:,f0idx(indf))/FWD_Cal_X(end,indf); %X_snaps = X_snaps/max(abs(X_snaps));
        Y_snaps = Y_p(:,f0idx(indf))/FWD_Cal_Y(end,indf); %Y_snaps = Y_snaps/max(abs(Y_snaps));
        Z_snaps = Z_p(:,f0idx(indf))/FWD_Cal_Y(end,indf); %Z_snaps = Z_snaps/max(abs(Z_snaps));
        Snaps   = [Hn_snaps;X_snaps;Y_snaps;Z_snaps];
        
        VS_BR = VS_W'*Snaps/384; VS_BR = reshape(VS_BR,l_phi,l_theta);
        VS_DEvsAZvsF(:,:,indf,block) = VS_BR.*conj(VS_BR);
        end
    end % end of frequency loop
        
    
    % read next chunk of data.
    VSData  = [VSData(:,OverLapPts +1: nPts)...
                fread(fidin,[Nchans OverLapPts],'float')];
    
end

%%
profile report
MAW = ones(1,10)/10;

thetaslice = find(theta==98);
fplot = find(f0==2000);

H_waterfall = filter(MAW,1,squeeze(H_DEvsAZvsF(:,thetaslice,fplot,:)),[],2);
H_waterfall = flipdim(H_waterfall,1);
VS_waterfall = filter(MAW,1,squeeze(VS_DEvsAZvsF(:,thetaslice,fplot,:)),[],2);
VS_waterfall = flipdim(VS_waterfall,1);

% plots for Azimuth slice response
H_FRAZ = squeeze(H_DEvsAZvsF(:,thetaslice,:,:));
VS_FRAZ = squeeze(VS_DEvsAZvsF(:,thetaslice,:,:));

figure(1)
subplot(211)
plot(phi,10*log10(abs(H_FRAZ(:,fplot))))
title(['H Only Azimuth response at ' num2str(f0(fplot)) ' Hz'],'FontSize',18)

subplot(212)
plot(phi,10*log10(abs(VS_FRAZ(:,fplot))))
title('VS Azimuth response','FontSize',18)
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('dB','FontSize',18)
xlim([phi(1) phi(end)])

%%
% plots for FRAZ responses
figure(2)
subplot(211)
imagesc(phi,f0,10*log10(abs(H_FRAZ(:,:,1))).')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Frequency, Hz','FontSize',18);
title(['H only FRAZ DE slice at $\theta$ = ' num2str(theta(thetaslice)) ' deg'],'FontSize',18)
colorbar
%     caxis([-40 0])

subplot(212)
imagesc(phi,f0,10*log10(abs(VS_FRAZ(:,:,1))).')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Frequency, Hz','FontSize',18);
title(['VS FRAZ DE slice at $\theta$ = ' num2str(theta(thetaslice)) ' deg'],'FontSize',18)
colorbar
%     caxis([-40 0])

%%
% plots for DE vs AZ responses
figure(3)
subplot(211)
H_DEvsAZ = squeeze(H_DEvsAZvsF(:,:,fplot,end));
imagesc(phi,theta,10*log10(abs(H_DEvsAZ)).')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('DE $\theta$, deg','FontSize',18);
title(['H only DE vs AZ at $f_o$ = ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
colorbar
%     caxis([-40 0])
subplot(212)
VS_DEvsAZ = squeeze(VS_DEvsAZvsF(:,:,fplot,end));
imagesc(phi,theta,10*log10(abs(VS_DEvsAZ)).')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('DE $\theta$, deg','FontSize',18);
title(['VS DE vs AZ at $f_o$ = ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
colorbar

drawnow

%%
figure(4)
subplot(211)
imagesc(phi,numBlocks:-1:1,H_waterfall.')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Time','FontSize',18);
title(['H only Waterfall at $f_o$ = ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
colorbar
subplot(212)
imagesc(phi,numBlocks:-1:1,VS_waterfall.')
colormap jet
xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Time','FontSize',18);
title(['VS Waterfall'],'FontSize',18)
colorbar