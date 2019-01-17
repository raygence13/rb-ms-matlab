% Make Single Signal for CSA example
close all
clear all
clc
present(0)

FINALVERSION = 1;
fs = 44.1e3;
f0 = 2e3;
c = 1500;
lambda = c/f0;
deltax = lambda/2;
Lsensors = 32;
% theta is radians from endfire
theta = 18*(pi/180); 
amplitude = 1;
Ns = ceil(50*(fs/f0)); % ten periods
s = sin((2*pi*f0/fs)*(0:Ns)).*tukeywin(Ns+1,1)';
s = [zeros(1,30) s zeros(1,70)];

figure(1)
plot((1:length(s))*(1/fs),s)
xlabel('Time (s)')
ylabel('Signal s[n]')

cleanx = MakeArrayData(s,deltax,Lsensors,c,fs,theta,amplitude);

SNR = -10;
sigmas2 = abs(s*s')/length(s);
sigman2 = sigmas2*10^(-SNR/10);
n = sqrt(sigman2)*randn(size(cleanx));
x = cleanx;

figure(2)
plot((1:size(x,2))*(1/fs),x);
xlabel('Time (s)')
ylabel('Signal x[n]')

if FINALVERSION
	save ArrayProject1CSA.mat x deltax Lsensors c fs
	save ArrayProject1CSAjbversion.mat
end

u = -1:0.01:1;
W = exp(-1i*(2*pi*f0/c)*deltax*(0:(Lsensors-1))'*u)./Lsensors;
Nfft = 2^nextpow2(length(x));
f = linspace(0,1-1/Nfft,Nfft).*fs;

X = fft(x,Nfft,2).*2/Nfft;
Xavg = sum(abs(X).^2,1)./sqrt(Lsensors);
myx = X(:,round((f0/fs)*Nfft));
Y = abs(W'*myx).^2;

figure(3)
plot(f./1e3,Xavg)
xlim([0 fs/2e3])
grid on
xlabel('f in kHz')

figure(4)
plot(u,Y)
