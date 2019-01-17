function [p_array, a]= GenVSETS(bbflag, AEL, phiS, thetaS)
% Generate vector sensor element time series
%
% generate simulated time series data from Wilcoxon vector sensor
%
% where:
%	bbflag: set to 1 to create broadband signal, else tone
%	AEL:    Array Element Locations [3 x #sensors]
%   phiS:   Source azimuth direction
%   thetaS: Source depth/elevation direction



if nargin<1,
    bbflag= 1;
end;
if size(AEL,1)~=3,
    AEL = AEL.';
end;
Nsensors = size(AEL,2);

% Environmental parameters
SL= 180;		% source level, dB referenced to 1 microPascal
% R= 100;			% range to source, meters
% Rho= 776;		% Isopar-M, kg per cubic meter
rho = 1026;
% C= 1280;		% Isopar-M, m/sec
c = 1500;

% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;


% Signal process parameters
bw= 500;		% broadband band width, Hz
upsample= 5;		% upsample factor for broadband
cf= 4000;		% center frequency of broadband
fs = 20e3;
%sr= upsample*bw;	% sample rate, Hz


T= 180;			% time period to generate, seconds

% generate pressure and acceleration signals
if bbflag,
    p0= randn(bw*T, 1);		% pressure wave, RMS value 1
    p0= p0*(10^(SL/20))/1000000;	% scale by source level
%     p0= p0/(R*R);			% scale by range, spherical spreading
    
    p1= resample(p0, 2*upsample, 1);	% upsample the broadband to 2x sample rate
    
    figure(1);
    specgram(p1, 2*upsample*bw, 1);
    
    p= p1.*cos(2*pi*cf*((0:length(p1)-1)')/(2*upsample*bw));	% shift to center frequency
    
    figure(2);
    specgram(p, 2*upsample*bw, 2);
    
    clear p0 p1;
    
else
    p= 1.414*cos(2*pi*cf*(1:(2*upsample*bw*T))/(2*upsample*bw))';	% generate tone
    
    p= p*(10^(SL/20))/1000000;	% scale by source level
%     p2= p2/(R*R);			% scale by range, spherical spreading
    
    figure(2);
    specgram(p, 2*upsample*bw, 2);
end;

uS = [cosd(phiS)*cosd(thetaS);
      sind(phiS)*sind(thetaS);
      sind(thetaS)];

tau = uS'*AEL/c;
Nfft = 2^nextpow2(length(p));

P = fft(p,Nfft);
P_Array = zeros(Nsensors,Nfft);

Tau(:,1:((Nfft/2))) = exp(-1i*(2*pi/Nfft*fs).*(tau'*(0:((Nfft/2)-1))));
Tau(:,((Nfft/2)+2):Nfft) = conj(flipdim(Tau(:,2:(Nfft/2)),2));


for sensor = 1:Nsensors
    P_Array(sensor,:) = Tau(sensor,:).*P.';
end
p_array = ifft(P_Array,Nfft,2);

% pmid= p_array(:,2:2:(end-1));					% for p, select mid-point each pair
pdot= (p_array(:,3:2:end)-p_array(:,1:2:(end-2)))*(upsample*bw);	% for pdot, calculate derivative
%   using sample each side of above
a= -pdot/(rho*c);					% convert pressure derivative to acceleration
%   in meters per second squared

% clear p2 pdot;
% 
% figure(3);
% specgram(pmid, upsample*bw, 3);
% 
% figure(4);
% specgram(a, upsample*bw, 4);
% 
% figure(5);
% subplot(2,1,1);
% plot(pmid(1:20), '-+');
% subplot(2,1,2);
% plot(a(1:20), '-+');
% 
% % model target position, azimuth and d/e
% de= 0;					% target in horizontal plane
% az= 2*pi*(0:(length(pmid)-1))'/length(pmid);	% circle around receiver, 0 to 2*pi (360 degrees)
% 
% ux= cos(az)*cos(de);
% uy= sin(az)*cos(de);
% uz= sin(de);
% 
% % observed physical parameters, m/sec and m/sec squared
% Ovs= [pmid, (ux.*a), (uy.*a), (uz.*a)];	% columns are vector sensor channels
% 
% figure(6);
% subplot(2,1,1);
% plot(Ovs(:,1));
% subplot(2,1,2);
% plot([Ovs(:,2), Ovs(:,3), Ovs(:,4)]);
% 
% % observed electrical vector sensor output with preAmp gain +20 dB, volts
% Vvs= 10*[1000000*(10^(Hs/20))*Ovs(:,1), Axs/9.8*Ovs(:,2), Ays/9.8*Ovs(:,3), Azs/9.8*Ovs(:,4)];
% 
% figure(7);
% plot(Vvs);
% 
% wavwrite(Vvs, upsample*bw, 'vs.wav');

% end;

