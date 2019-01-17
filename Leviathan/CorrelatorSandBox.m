% Leviathan Simulation  - R. Bautista 12 Feb 2018
% Generate time series data to evaluate performance of correlator detector

clear
close all

%% Define sensor locations in [m]
% pLX = [3 x 4] position vector for Leviathan XLUUV (square geometry)
% pSVD = [3 x 4] position vector for SVD (T geometry)
doplot = 0;
pLX = [-0.3048, 0.3048,  0;...
        0.3048, 0.3084,  0;...
        0.3048, -0.3084, 0;...
       -0.3048, -0.3084, 0]';
   
% pSVD = [-0.4572, 0.6096,  0;...
%          0.4572, 0.6096,  0;...
%          0,      0,       0;...
%          0,      -2.4384, 0]';

MaxDistance = sqrt(sum((pLX(:,1) - pLX(:,3)).^2));
if doplot
    figure
    subplot(121)
    plot(pLX(1,:),pLX(2,:),'.')
    subplot(122)
    plot(pSVD(1,:),pSVD(2,:),'.')
end

%% Generate time series of received signal
c = 1500;                   % sound speed [m/s]
f0 = floor(2*c/2.4384);     % frequency of signal based on T geometry spacing
fs = f0*10;                 % sampling frequency [Hz]
nSamples = 1000;

% Define signal direction
phiS = 85;      % [deg]
thetaS = 10;    % [deg]
uS = [sind(phiS)*cosd(thetaS);
      sind(phiS)*sind(thetaS);
      cosd(thetaS)];
tauLX = uS'*pLX/c;   % time delay [s]
% tauSVD = uS'*pSVD/c; % time delay [s]

Nperiods = 50;
A = 1;  % Amplitude
Ns = ceil(Nperiods*(fs/f0)); % number of samples required to generate Nperiods
s = A*cos((2*pi*f0/fs)*(0:Ns)).*tukeywin(Ns+1,1)';
s = [zeros(1,round(MaxDistance*20/c*fs)) s zeros(1,round(MaxDistance*20/c*fs))];
Nfft = 2^nextpow2(length(s));

% plot(0:length(s)-1,s)


S = fft(s,Nfft);
S_LX = zeros(4,Nfft);
% S_SVD = zeros(4,Nfft);

Tau_LX(:,1:((Nfft/2))) = exp(-1i*(2*pi/Nfft*fs).*(tauLX'*(0:((Nfft/2)-1))));
Tau_LX(:,((Nfft/2)+2):Nfft) = conj(flipdim(Tau_LX(:,2:(Nfft/2)),2));

% Tau_SVD(:,1:((Nfft/2))) = exp(-1i*(2*pi/Nfft*fs).*(tauSVD'*(0:((Nfft/2)-1))));
% Tau_SVD(:,((Nfft/2)+2):Nfft) = conj(flipdim(Tau_LX(:,2:(Nfft/2)),2));

for sensor = 1:4
    S_LX(sensor,:) = Tau_LX(sensor,:).*S;
    
%     S_SVD(sensor,:) = Tau_SVD(sensor,:).*S;
end
s_LX = ifft(S_LX,Nfft,2);
% s_SVD = ifft(S_SVD,Nfft,2);

SNR = 0;   % Desired SNR [dB]
varW = A/10^(SNR/10)/2/Nfft;
n = varW.*randn(4,Nfft,nSamples);
r_ts_LX = bsxfun(@plus,s_LX,n);
% r_ts_SVD = bsxfun(@plus,s_SVD,n);
figure
plot((0:length(r_ts_LX)-1)/fs,r_ts_LX(:,:,1))

R_f_LX = fft(r_ts_LX,Nfft,2);
% R_f_SVD = fft(r_ts_SVD,Nfft,2);
%% Generate signal with A~CN(0,varS)
varS = 1;
As = sqrt(varS/2).*(randn(1,nSamples) + 1i.*randn(1,nSamples));
v_LX = exp(-1i*2*pi*f0*tauLX);
snaps_LX = v_LX'*As;

% v_SVD = exp(-1i*2*pi*f0*tauSVD);
% snaps_SVD = v_SVD'*As;

varW = 1;
snaps_n = sqrt(varW/2).*(randn(4,nSamples) + 1i.*randn(4,nSamples));
