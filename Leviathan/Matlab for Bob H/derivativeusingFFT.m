function [t] = derivativeusingFFT(N, Signal)
% derivativeusingFFT(N, Signal)
%
% where:  N       is number of samples of waveform for demo
%         Signal  1-> tone, 2-> 3 tones, 3-> noise
%

Fs= 25000;                  % sample rate 25 kHz
f= 1002;                    % tone at 1 kHz

t= (0:(N-1))/Fs;  % time at sample points

% gain portion to apply to each frequency point
g0= 0:(N/2);
g1= -(((N/2)-1):-1:1);
g= [g0 g1] * 2*pi*Fs/N;
figure(1);
plot([g1 g0]);

% sample wave, tone
if (Signal==1),
 s= sin(2*pi*f*t);
end;

% sample wave, 3 tones
if (Signal==2),
 s= sin(2*pi*f*t)+.8*sin(2*pi*(f+50)*t)-1.2*sin(2*pi*(f+75)*t);
end;

% sample wave, filtered noise
if (Signal==3),
 FO= 30;
 Ctr= .667;
 sz=rand(1,N+4*FO)-.5;
 fir=firls(30, [0 Ctr-.05 Ctr+.05 1], [1 1 0 0]);
 s=filter(fir, 1, sz);
 s=s((2*FO):(N+2*30-1));
end;


h= [0 hanning(N-2)' 0];
h= [0*(1:10) 1+0*(1:(N-20)) 0*(1:10)];

h0= .5*(1-cos(pi*(0:9)/10));
h1= .5*(1+cos(pi*(1:10)/10));
h= [h0 1+0*(1:(N-20)) h1];

s= s.*h;

% convert to frequency domain
S= fft(s);

% take derivative in frequency domain, apply gain
SD= S.*g;

% then apply 90 degree phase shift (switch real/imag parts)
SD= i*real(SD)-imag(SD);

% convert back to time domain
sd= ifft(SD);
sdi= 0:(N-1);

% take derivative using difference sample to sample
sd0= diff(s)*Fs;
sd0i= (0:(N-2))+.5;

sd1= (s(2:N)-s(1:(N-1)))*Fs;

figure(2);
plot(sdi, [real(sd)' imag(sd)'], sd0i, sd0, '+', sd0i, sd1, 'o');

end
