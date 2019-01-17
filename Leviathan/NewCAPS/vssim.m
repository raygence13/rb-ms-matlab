function ts= vssim(bbflag, num_sensors)
% ts= vssim(bbflag, num_sensors)
%
% generate simulated time series data from Wilcoxon vector sensor
%
% where:
%	bbflag			set to 1 to create broadband signal, else tone
%	num_sensors		no. of sensors to simulate, 4 channels per sensor
%
%

close all;
ts= [];

if nargin<1,
 bbflag= 1;
 end;
if nargin<2,
 num_sensors= 1;
 end;

% Environmental parameters
SL= 180;		% source level, dB referenced to 1 microPascal
R= 100;			% range to source, meters
Rho= 776;		% Isopar-M, kg per cubic meter
C= 1280;		% Isopar-M, m/sec


% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;


% Signal process parameters
bw= 5000;		% broadband band width, Hz
upsample= 5;		% upsample factor for broadband
cf= 4000;		% center frequency of broadband
%sr= upsample*bw;	% sample rate, Hz


T= 180;			% time period to generate, seconds

% generate pressure and acceleration signals
if bbflag,
 p0= randn(bw*T, 1);		% pressure wave, RMS value 1
 p0= p0*(10^(SL/20))/1000000;	% scale by source level
 p0= p0/(R*R);			% scale by range, spherical spreading

 p1= resample(p0, 2*upsample, 1);	% upsample the broadband to 2x sample rate

 figure(1);
 specgram(p1, 2*upsample*bw, 1);

 p2= p1.*cos(2*pi*cf*((1:length(p1))')/(2*upsample*bw));	% shift to center frequency

 figure(2);
 specgram(p2, 2*upsample*bw, 2);

 clear p0 p1;

else,
 p2= 1.414*cos(2*pi*cf*(1:(2*upsample*bw*T))/(2*upsample*bw))';	% generate tone

 p2= p2*(10^(SL/20))/1000000;	% scale by source level
 p2= p2/(R*R);			% scale by range, spherical spreading

 figure(2);
 specgram(p2, 2*upsample*bw, 2);
 end;

p= p2(2:2:(end-1));					% for p, select mid-point each pair
pdot= (p2(3:2:end)-p2(1:2:(end-2)))*(upsample*bw);	% for pdot, calculate derivative
							%   using sample each side of above
a= -pdot/(Rho*C);					% convert pressure derivative to acceleration
							%   in meters per second squared

clear p2 pdot;

figure(3);
specgram(p, upsample*bw, 3);

figure(4);
specgram(a, upsample*bw, 4);

figure(5);
subplot(2,1,1);
plot(p(1:20), '-+');
subplot(2,1,2);
plot(a(1:20), '-+');

% model target position, azimuth and d/e
de= 0;					% target in horizontal plane
az= 2*pi*(0:(length(p)-1))'/length(p);	% circle around receiver, 0 to 2*pi (360 degrees)

ux= cos(az)*cos(de);
uy= sin(az)*cos(de);
uz= sin(de);

% observed physical parameters, m/sec and m/sec squared
Ovs= [p, (ux.*a), (uy.*a), (uz.*a)];	% columns are vector sensor channels

figure(6);
subplot(2,1,1);
plot(Ovs(:,1));
subplot(2,1,2);
plot([Ovs(:,2), Ovs(:,3), Ovs(:,4)]);

% observed electrical vector sensor output with preAmp gain +20 dB, volts
Vvs= 10*[1000000*(10^(Hs/20))*Ovs(:,1), Axs/9.8*Ovs(:,2), Ays/9.8*Ovs(:,3), Azs/9.8*Ovs(:,4)];

figure(7);
plot(Vvs);

wavwrite(Vvs, upsample*bw, 'vs.wav');


