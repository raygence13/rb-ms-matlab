function [wf, iwf]=vscbf(baseband, wavfile)
% [wf, iwf]=vscbf(baseband, wavfile)
%
% read wavfile of vector sensor time series data and form beams
%
% where:
%	baseband -> 1, translate carrier (center frequency) to DC
%	wavfile	filename of .wav file
%

close all;
wf= [];
iwf= [];

if nargin<1,
 baseband= 0;
 end;
if nargin<2,
 wavfile= 'vs.wav';
 end;


% Signal process parameters
Tavg0= 0.5;	% desired averaging interval

% Environmental parameters
Rho= 776;		% Isopar-M, kg per cubic meter
C= 1280;		% Isopar-M, m/sec

% Sensor parameters
Hs= -165;		% dB // 1 Volt per microPascal
Axs= 1.5;		% Volts per G
Ays= 1.5;
Azs= 1.5;


[Vvs, sr0]= wavread(wavfile);

% Signal parameters
bw= 5000;		% broadband band width, Hz
upsample= 5;		% upsample factor for broadband
cf= 4000;		% center frequency of broadband


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
pause(0);


if baseband,
 % compute carrier modulation to shift positive frequencies at cf to 0 (DC)
 carrier= exp(-j*2*pi*cf*((1:length(H))')/sr0);

 % apply modulation to real signal, positive frequencies centered on 0 Hz,
 % negative frequencies shift further negative, perhaps wrapping around to positive
 H= H.*carrier;
 Ax= Ax.*carrier;
 Ay= Ay.*carrier;
 Az= Az.*carrier;

 % downsample by 4 and lowpass filter at cf/4 Hz, result is complex waveform
 H= resample(H, 1, 4);
 Ax= resample(Ax, 1, 4);
 Ay= resample(Ay, 1, 4);
 Az= resample(Az, 1, 4);

 % adjust sample rate
 sr= sr0/4;

else,
 sr= sr0;
 end;

figure(3);
specgram(H, sr, 3);
pause(0);


% fft size in seconds, integer samples nearest power of 2 larger
fftlen= 2^ceil(log2(.1*sr));

% fft size in time, frequency of FFT bins
Tfft= fftlen/sr;	% duration of integer chunk
df= 1/Tfft;		% frequency increment of FFT

% be sure averaging time (Tavg0) selected above is greater than FFT time
if Tavg0<Tfft,
 Tavg0= Tfft;
 end;

if baseband,
 disp(sprintf('Processing complex data with band center %4d shifted to DC', cf));
else,
 disp(sprintf('Processing real data with band center at %4d', cf));
 end;
disp(sprintf('FFT size: %5d samples, %6.3f sec; Power average %6.3f sec', fftlen, Tfft, Tavg0));

% derivative exp(jwt) in frequency domain multiplies by jw*exp(jwt)
jw= [df*((1:(fftlen/2))-1), 0, -df*((fftlen/2)-(2:(fftlen/2))+1)]';
if baseband,
 jw= j*2*pi*(jw+cf);
else,
 jw= j*2*pi*(jw+0);
 end;


% averaging interval, no. of averages and time period
navg= floor(Tavg0/Tfft);
if navg<2,
 navg= 2;
 end;
Tavg= navg*Tfft;
disp(sprintf('Average time: desired %5.2f, actual %5.2f seconds, %2d counts', Tavg0, Tavg, navg));

% calculate steering vectors, 0 to 360 in .5 degree increments
v= [];
for bm=0:0.5:360,
 az= bm;
 de= 0;
 v= [v, [ 1; cos(pi*az/180)*cos(pi*de/180); sin(pi*az/180)*cos(pi*de/180); sin(pi*de/180)]];
 end;

fftstart= 0;

avg= 0;
bmaccum= 0*(0:0.5:360);
ibmaccum= 0*(0:0.5:360);

% beamforming loop, do chunks
while (fftstart+fftlen)<=length(H),
 % selection index
 fftindex= fftstart+(1:fftlen);

 % select time segment of sensor data
 h= H(fftindex);
 x= Ax(fftindex);
 y= Ay(fftindex);
 z= Az(fftindex);

 % convert hydrophone pressure/velocity to acceleration
 hdot= ifft(jw.*fft(h));	% derivative in frequency domain
 hdot= -hdot/(Rho*C);		% scale by Rho*C
 %hdot= real(hdot);		% take real part

 %figure(3);
 %plot([hdot, x, y, z]);

 % SSM
 ssm= ([hdot, x, y, z]')*[hdot, x, y, z];

 % scale by number of time samples
 ssm= ssm./size(hdot, 1);

 % loop over beams, conventional beamforming
 for bm=1:size(v, 2),
  bmaccum(bm)+= (v(:,bm)')*ssm*v(:,bm);
  end;

 % take inverse of SSM above
 [issm, r]= inv(ssm(1:3,1:3));

 % check inverse OK, original not singular
 if r>0,
  % loop over beams, apply MPDR Adaptive beamforming
  for bm=1:size(v, 2),
   ibmaccum(bm)+= 1/((v(1:3,bm)')*issm*v(1:3,bm));
   end;
  end;

 avg= avg+1;

 if (avg==navg),
  bmaccum= bmaccum/navg;
  ibmaccum= ibmaccum/navg;

  figure(4);
  subplot(2,1,1);
  plot(real(bmaccum));
  subplot(2,1,2);
  plot(real(ibmaccum));

  wf= [wf; real(bmaccum)];
  iwf= [iwf; real(ibmaccum)];

  avg= 0;
  bmaccum= 0*bmaccum;
  ibmaccum= 0*ibmaccum;
  end;


 fftstart= fftstart+fftlen;
 pause(0);
 endwhile;

figure(5);
imagesc(10*log10(wf));
figure(6);
imagesc(10*log10(iwf));

end;

