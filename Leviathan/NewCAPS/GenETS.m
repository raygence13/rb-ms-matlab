function [ETS, p] = GenETS(bbflag,depth,r,phiS)
%Function to generate clean signal element time series for hydrophones and
%vector sensors.
% RB 9MAR2018. RB 15MAR2018.
% Inputs:
% bbflag:   option to generate broadband signal
% depth:    array depth (negative since +z is up) [m]
% r:        slant range of target
% phiS:     target azimuth bearing (deg)
%
% Outputs:
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

% check to see if range is greater than depth
if r < abs(depth)
    r = abs(depth) + 0.1*abs(depth);
    warning(['Range must be greater than depth. Range changed to ' num2str(r)])
end

% define array element locations [3 x p.Nsensors]
c   = 1500;       % sound speed seawater [m/s]
rho = 1026;     % seawater density [kg/m^3]
p   = ArrayElementLocations(c,rho,depth);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define signal time series parameters
fs          = p.f0*20; % sampling frequency [Hz]
Nperiods    = 100;   % number of periods desired in single tone signal
Ns          = ceil(Nperiods*(fs/p.f0)); % number of samples required to generate Nperiods for p.f0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~bbflag
    % Create single tone signal
    sig = cos((2*pi*p.f0/fs)*(0:Ns-1));    % generate single tone
else
    % Create bandlimited signal
    BW = 5e3;   % bandwidth [Hz]
    B   = fir1(100,(BW-50)/fs,'low');
    B   = B/sum(B);
    fgwn = filter(B,1,randn(1,Ns));    % filtered Gaussian white noise
    sig = fgwn.*cos(2*pi*2550/fs.*(0:Ns-1));    % shift in frequency
    clear D Filter fgwn
end

% Define target parameters
SL      = 180;   % source level, dB referenced to 1 microPascal
thetaS  = -round(asind(p.AEL(3,3)/r)*10)/10;   % rounded target de bearing based on nominal range [deg]
A       = sqrt(2)*(10^(SL/20)*1e-6)/r^2;  % Amplitude with RMS = 1, with range attenuation and source level

% target unit vector
uS = [cosd(phiS)*cosd(thetaS);
      sind(phiS)*cosd(thetaS);
      sind(thetaS)];
% time delays based on target location and p.AEL [1 x p.Nsensors]
tau = uS'*p.AEL/c;

sig = A*[zeros(1,round(sum(abs(tau))/2*fs)) sig zeros(1,round(sum(abs(tau))/2*fs))];
% pad signal with zeros to allow for time delays
Nfft = 2^nextpow2(length(sig));                 % fft size
frequencies = linspace(0,1-1/Nfft,Nfft).*fs;    % frequency vector
Sig = fft(sig,Nfft);                            % signal spectrum


% Initialize memory for hydrophone data [Nfft x p.Nsensors]
H = zeros(Nfft,p.Nsensors);

% time delays in freq domain [p.Nsensors x Nfft]
Tau(:,1:((Nfft/2)))       = exp(-1i*(2*pi/Nfft*fs).*(tau'*(0:((Nfft/2)-1))));
Tau(:,((Nfft/2)+2):Nfft)  = conj(flipdim(Tau(:,2:(Nfft/2)),2));

% applying time delays to signal in freq domain
for sensor = 1:p.Nsensors
    H(:,sensor) = Tau(sensor,:).*Sig;
end

% time series hydrophone data
hData     = ifft(H,Nfft,1);
time      = (0:length(hData)-1)/fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(size(H)); % Initialize memory for acceleration spectrum
% Convert hydrophone pressure time series to velocity by dividing by rho*c
% Convert hydrophone velocity to acceleration via derivative in freq. domain by multiplying by jw
A(1:Nfft/2,:)         = bsxfun(@times,H(1:Nfft/2,:)./(rho*c),1i*2*pi.*frequencies(1:Nfft/2).');
A((Nfft/2)+2:Nfft,:)  = bsxfun(@times,H(Nfft/2+2:Nfft,:)./(rho*c),-1i*2*pi*frequencies(2:Nfft/2).');

a     = ifft(A,Nfft,1); % acceleration time series
xData = a*uS(1);
yData = a*uS(2);
zData = a*uS(3);

% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;

% observed electrical vector sensor output with preAmp gain +20 dB, volts
ETS.h       = 10*(10^(Hs/20)/1e-6)*hData;
ETS.Hs      = Hs;
ETS.a       = a;
ETS.x       = 10*Axs/9.8*xData;
ETS.Axs     = Axs;
ETS.y       = 10*Ays/9.8*yData;
ETS.Ays     = Ays;
ETS.z       = 10*Azs/9.8*zData;
ETS.Azs     = Azs;
ETS.time    = time;
ETS.fs      = fs;
ETS.Nfft    = Nfft;
ETS.frequencies = frequencies;
ETS.f0      = p.f0;
ETS.az      = phiS;
ETS.de      = thetaS;
ETS.range   = r;

% figure('position',[1950 10 1600 800])
% subplot(211)
% plot(ETS.time,ETS.h)
% xlim([ETS.time(1) ETS.time(end)]); xlabel('t [s]')
% if ~bbflag
%     title(['Narrow band signal time series, $f_0 = $' num2str(ETS.f0) ' Hz'])
% else
%     title(['Broad band signal time series, with ' num2str(BW) ' Hz BW'])
% end
% subplot(212)
% plot(ETS.frequencies(1:Nfft/2),abs(Sig(1:Nfft/2)))
% xlim([ETS.frequencies(1) ETS.frequencies(Nfft/2)]);
% xlabel('f [Hz]')
% title(['Signal Spectrum'])
end


function p = ArrayElementLocations(c,rho,depth)
% Cartesian coordinate system:
% +x: forward,      phiS = 0,  theta = 0,       ux = cos(phi)cos(theta)
% +y: starboard,    phiS = 90, theta = 0,       uy = sin(phi)cos(theta)
% +z: down,                    thetaS = +90,    uz = sin(theta)

% define array element locations [3 x Nsensors] in units of [m]
% geometry defined to be square
% depth = 300;        % depth of array, meters
dx = (45/12)*0.3048;       % sensor spacing in [m] (ft x .3048m/ft)
dy = (43/12)*0.3048;
AEL = [0, 0, depth;...
       0, -dy, depth;...
       -dx, 0, depth;...
       -dx, -dy, depth].'; 
% AEL = [d*(0:9)', zeros(10,1), -depth*ones(10,1)]';    % line array geometry to double check VS response

% center array around x-y origin
% AEL(1,:) = AEL(1,:) - mean(mean(AEL(1,:)));
% AEL(2,:) = AEL(2,:) - mean(mean(AEL(2,:)));
Nsensors = size(AEL,2);     % number of sensors
MaxDistance = sqrt(sum((AEL(:,1) - AEL(:,4)).^2));  % find max distance between pair of sensors [m]
lambda0 = dy*2;
f0 = floor((c/lambda0)/10)*10;

p.AEL = AEL;
p.Nsensors = Nsensors;
p.MaxDistance = MaxDistance;
p.depth = depth;
p.c = c;
p.rho = rho;
p.f0 = f0;

% figure('position',[1950 10 1600 800])
% plot(p.AEL(2,:)/0.3048,p.AEL(1,:)/0.3048,'o','MarkerSize',20)
% xlim([-dx/0.3048-0.1 0+0.1]); ylim([-dx/0.3048-0.1 0+0.1])
% axis equal
% % title('BASS Array Element Locations')
% ylabel({'[ft]','AFT $\leftarrow$y$\rightarrow$ FWD'}); xlabel({'PORT $\leftarrow$x$\rightarrow$ STBD','[ft]'});
end
