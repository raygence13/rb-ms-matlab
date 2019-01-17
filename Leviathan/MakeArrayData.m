function [ x ] = MakeArrayData(p,M,N,f0,c,fs,theta,phi,amplitude)
%MAKEARRAYDATA created array time series data
%   Creates time series for each sensor in the array for the arriving
%   signal s for a uniform linear array of Nsensors separated by deltax.
%   The propagation speed is c, the sampling frequency is fs.  The sources
%   are at direction of arrival given in theta (radians) relative to
%   endfire, so 0 to pi.  Each arrival is scaled by the amplitude given in
%   that vector (default = 1).  If amplitude is a scalar, all of the
%   signals have the same amplitude given by amplitude.   
%   The output x is Nsensors by Ntime, so each row is the time series for
%   one sensor and Ntime exceeds the length of s by the amount of time
%   needed to make a power of 2 above the maximum array transit 

% john buck
% 14 july 2014

Nperiods = 50;
Ns = ceil(Nperiods*(fs/f0)); % ten periods
s = sin((2*pi*f0/fs)*(0:Ns)).*tukeywin(Ns+1,1)';
s = [zeros(1,30) s zeros(1,70)];

Ns = size(s,2);
% propagation time across the array
MaxDelay = abs(p(2)-p(end-1))/c*fs*100;
Nfft = 2^nextpow2(Ns+MaxDelay);

u = [sin(phi)*cos(theta);
    sin(phi)*sin(theta);
    cos(theta)];

taun = (u.'*p)./c*fs;
x = zeros(M*N,Nfft);

% only one signal
if (size(s,1)==1)
	S = fft(s,Nfft);
	W = zeros(1,Nfft);
	% note leave W(Nfft/2)+1) as zero since this is omega = pi;
	
	for xidx = 1:M*N
        W(1:((Nfft/2))) = exp(-1i*(2*pi/Nfft)*taun(xidx)*(0:((Nfft/2)-1)));
        W(((Nfft/2)+2):Nfft) = conj(fliplr(W(2:(Nfft/2))));
		temp = S.*W;
		x(xidx,:) = ifft(temp);
	end
	x = x*amplitude(1);
else
	error(mfilename,'Lazy jb needs to write this case still');
end

		
	


return

