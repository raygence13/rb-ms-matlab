function ts=vsibf(wavfile)
% ts=vsibf(wavfile)
%
% read wavfile of vector sensor time series data and form beams
%
% where:
%	wavfile	filename of .wav file
%

close all;
ts= [];

if nargin<1,
 wavfile= 'vs.wav';
 end;


% Environmental parameters
Rho= 776;		% Isopar-M, kg per cubic meter
C= 1280;		% Isopar-M, m/sec

% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;


[Vvs, sr]= wavread(wavfile);

% backout preAmp gain
Vvs= Vvs/10;

figure(1);
plot(Vvs);

% backout sensitivities
H= Vvs(:,1)/(10^(Hs/20)*1000000);	% pressure in Pascals
Ax= Vvs(:,2)*9.8/Axs;			% X component of acceleration, meter per second squared
Ay= Vvs(:,3)*9.8/Ays;			% Y component of acceleration, meter per second squared
Az= Vvs(:,4)*9.8/Azs;			% Z component of acceleration, meter per second squared

figure(2);
subplot(2,1,1);
plot(H);
subplot(2,1,2);
plot([Ax, Ay, Az]);

% chunk size (averaging time) in seconds, integer samples
T= .5;
chunklen= ceil(sr*T);

% be sure even
chunklen= 2*floor(chunklen/2);


% actual chunk length in time
T0= chunklen/sr;	% duration of integer chunk
df= 1/T0;		% frequency increment of FFT

% derivative exp(jwt) in frequency domain multiplies by jw*exp(jwt)
jw= [2*pi*df*((1:(chunklen/2))-1), 0, -2*pi*df*((chunklen/2)-(2:(chunklen/2))+1)]';
jw= j*jw;



chunkstart= 0;
brg= [];

% beamforming loop, do chunks
while (chunkstart+chunklen)<=length(H),
 % selection index
 chunkindex= chunkstart+(1:chunklen);

 % select time segment of sensor data
 h= H(chunkindex);
 x= Ax(chunkindex);
 y= Ay(chunkindex);
 z= Az(chunkindex);

 
 % convert hydrophone pressure/velocity to acceleration
 hdot= ifft(jw.*fft(h));	% derivative in frequency domain
 hdot= -hdot/(Rho*C);		% scale by Rho*C
 hdot= real(hdot);		% take real part

 figure(3);
 plot([hdot, x, y, z]);

 % calculate direction vector using intensity, good only for 1 target and minimal noise
 s= ([x, y, z]')*hdot;
 u= s/sqrt(sum(s.*s));

 % calculate azimuth bearing and d/e
 az= atan2(u(2), u(1))*180/pi;
 az= mod(az, 360);
 de= asin(u(3))*180/pi;
 brg= [brg, az];

 disp(sprintf('Az %8.1f, D/ %6.1f', az, de));
  
 pause(0);

 chunkstart= chunkstart+chunklen;
 endwhile;

figure(4);
plot(brg, '+');

end;

