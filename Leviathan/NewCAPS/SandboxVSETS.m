clear
close all
% define array element locations [3 x Nsensors]
AEL = [-0.3048, 0.3048,  0;...
        0.3048, 0.3084,  0;...
        0.3048, -0.3084, 0;...
       -0.3048, -0.3084, 0]';
% find max distance between pair of sensors
MaxDistance = sqrt(sum((AEL(:,1) - AEL(:,3)).^2));
Nsensors = size(AEL,2);     % number of sensors

c = 1500;       % sound speed seawater [m/s]
rho = 1026;     % seawater density [kg/m^3]
% Cartesian coordinate system:
% +x: forward, phiS = 0
% +y: starboard
% +z: down, thetaS = +90

profile on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Signal process parameters (some from B. Hayford vssim code)
SL = 180;           % source level, dB referenced to 1 microPascal
r = 100;			% range to source, meters
% SampleSize = 100;
% bw = 1000;		% broadband bandwidth [Hz]
% T = 1;         % time period to generate [sec]
% upsample = 5;	% upsample factor for broadband
% f0 = 5000;		% center frequency of broadband [Hz]
% fs = 2*upsample*(bw+f0);      % sampling rate [Hz]
% 
% p0 = randn((bw+f0)*T,SampleSize)*(10^(SL/20)*1e-6)/r^2;   % pressure wave, RMS value 1, scaled by source level, scaled  by range
% p1 = resample(p0,2*upsample,1);	% upsample the broadband to 2x sample rate
% signal = bsxfun(@times,p1,cos(2*pi*f0*((0:length(p1)-1)'/fs)));	% shift to center frequency
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create single tone signal
f0 = floor(2*c/MaxDistance);    % frequency of signal based on largest spacing
fs = f0*10;                     % sampling frequency [Hz]
Nperiods = 100;
Ns = ceil(Nperiods*(fs/f0)); % number of samples required to generate Nperiods
A = sqrt(2)*(10^(SL/20)*1e-6)/r^2;  % Amplitude with range attenuation and source level

sig = A*cos((2*pi*f0/fs)*(0:Ns)).*tukeywin(Ns+1,1)';    % generate single tone
sig = [zeros(1,round(MaxDistance*20/c*fs)) sig zeros(1,round(MaxDistance*20/c*fs))];
% pad signal with zeros to allow for time delays
Nfft = 2^nextpow2(length(sig)); % fft size
Sig = fft(sig,Nfft);    % signal spectrum

phiS = 30;       % target az bearing [deg]
thetaS = -30;   % target de bearing [deg]
% target unit vector
uS = [cosd(phiS)*cosd(thetaS);
      sind(phiS)*cosd(thetaS);
      sind(thetaS)];
% time delays based on target location and AEL [1 x Nsensors]
tau = uS'*AEL/c;

% Initialize memory for hydrophone data [Nfft x Nsensors]
H = zeros(Nfft,Nsensors);
% time delays in freq domain [Nsensors x Nfft]
Tau(:,1:((Nfft/2))) = exp(-1i*(2*pi/Nfft*fs).*(tau'*(0:((Nfft/2)-1))));
Tau(:,((Nfft/2)+2):Nfft) = conj(flipdim(Tau(:,2:(Nfft/2)),2));
% applying time delays to signal in freq domain
for sensor = 1:Nsensors
    H(:,sensor) = Tau(sensor,:).*Sig;
end
% time series hydrophone data
h = ifft(H,Nfft,1);
time = (0:length(h)-1)/fs;

% add noise to hydrophone signal
varP = trace(h*h')/Nsensors;
SNR = 20;
varNp = varP/(10^(SNR/10));

SampleSize = 1000;
hnoise = varNp.*randn(Nfft,Nsensors,SampleSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. Hayford code from vssim to generate accelerometer time series
dh = h(2:end,:)-h(1:end-1,:);
dt = time(2:end)-time(1:end-1);
hdot= bsxfun(@times,dh,1./dt.');	% for pdot, calculate derivative
 % using sample each side of above
a= -hdot/(rho*c);
varA = trace(a*a')/Nsensors;
varNa = varA/(10^(SNR/10));
anoise = varA.*randn(length(a),Nsensors,SampleSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hData = bsxfun(@plus,h,hnoise);
aData = bsxfun(@plus,a,anoise);

% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;

xData = aData*uS(1);
yData = aData*uS(2);
zData = aData*uS(3);

% observed electrical vector sensor output with preAmp gain +20 dB, volts
hData = 10*(10^(Hs/20)/1e-6)*hData;
xData = 10*Axs/9.8*xData;
yData = 10*Axs/9.8*yData;
zData = 10*Axs/9.8*zData;

figure(1)
plot(time,hData(:,1,1),(0:length(a)-1)/fs,aData(:,1,1))
profile report


%% Beamforming

% 3-dimensional beamforming parameters
phi = -180:180;      % degrees
l_phi = length(phi);
theta = -90:0;
l_theta = length(theta);

ux = cosd(phi)'*cosd(theta);   % sensors along z-axis
uy = sind(phi)'*cosd(theta);
uz = ones(l_phi,1)*sind(theta);

ux = ux(:);     % x-component of unit vector
uy = uy(:);     % z-component of unit vector
uz = uz(:);     % y-component of unit vector

% 3-dimensional unit vector
u = zeros(3,l_phi*l_theta);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);


HData = fft(hData,Nfft,1);
XData = fft(xData,Nfft,1);
YData = fft(yData,Nfft,1);
ZData = fft(zData,Nfft,1);

