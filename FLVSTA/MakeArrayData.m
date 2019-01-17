function [ h_data] = MakeArrayData(p,Nchans,f0,c,fs,theta_tgt,phi_tgt,amplitude)
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

Ns = ceil(100*(fs/f0)); % 100 periods of signal
s = sin((2*pi*f0/fs)*(0:Ns)).*tukeywin(Ns+1,1)';

% propagation time across the array
MaxDelay = round(abs(p(1)-p(end))/c*fs*50);
s = [zeros(1,MaxDelay) s zeros(1,MaxDelay)];
Nfft = 2^nextpow2(length(s));

u = [cos(theta_tgt);
    sin(theta_tgt)*sin(phi_tgt);
    sin(theta_tgt)*cos(phi_tgt)];

% u = [sin(phi_tgt)*cos(theta_tgt);
%     cos(phi_tgt)*cos(theta_tgt);
%     sin(theta_tgt)];

taun = -(u.'*p)./c*fs;
pressure_data = zeros(Nchans,Nfft);

% only one signal
if (size(s,1)==1)
	S = fft(s,Nfft);
	W = zeros(Nchans,Nfft);
    W(:,1:((Nfft/2))) = exp(-1i*(2*pi/Nfft)*taun.'*(0:((Nfft/2)-1)));
    W(:,((Nfft/2)+2):Nfft) = conj(fliplr(W(:,2:(Nfft/2))));
	% note leave W(Nfft/2)+1) as zero since this is omega = pi;
	
	for xidx = 1:Nchans
		temp = S.*W(xidx,:);
		pressure_data(xidx,:) = ifft(temp,Nfft);
	end
	h_data = pressure_data*amplitude(1);
%     x_data = pressure_data*u;
%     y_data = pressure_data*u(2);
%     z_data = pressure_data*u(3);
    
else
	error(mfilename,'Lazy jb needs to write this case still');
end

		
	


return

